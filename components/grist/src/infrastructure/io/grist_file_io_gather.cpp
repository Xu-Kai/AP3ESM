#include <mpi.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <stdexcept>
#include <pnetcdf.h>
#include <algorithm>
#include <regex>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

extern "C" void create_directory(const char *parten_dir, const char *filename, const char *varname,
                                 char *outpath, int outlen)
{
    // 1. 定义正则
    regex pattern(R"((\d{4}-\d{2}-\d{2}-\d+)\.(.+)\.nc)");
    cmatch match;
    string final_path;

    // printf("Input filename: %s\n", filename);
    // printf("Creating directory for variable: %s\n", varname);
    // printf("Input parent dir: %s\n", parten_dir);
    // cout << "Creating directory for variable: " << varname << endl;
    // cout << "Parent dir: " << parten_dir << ", filename: " << filename << endl;

    // 2. 查找匹配
    if (regex_search(filename, match, pattern))
    {
        string date_str = match[1].str();
        string suffix_str = match[2].str();

        string path1 = parten_dir + date_str;
        mkdir(path1.c_str(), 0777);

        string path2 = path1 + "/" + string(varname) + ".dir";
        mkdir(path2.c_str(), 0777);

        // string final_path = path2 + "/" + suffix_str;
        final_path = path1 + "/";
    }
    else
    {
        final_path.clear();
    }

    int n = static_cast<int>(final_path.size());
    if (n > outlen - 1)
        n = outlen - 1;
    if (n > 0)
        memcpy(outpath, final_path.data(), static_cast<size_t>(n));
    if (outlen > 0)
        outpath[n] = '\0';
}

extern "C" void wrap_output_1d_gather_sp_c(const float *temp, int local_len,
                                           int color_id, int all_comm_f,
                                           const char *outdir, const char *filename,
                                           int open_mode,
                                           int my_count, int my_start,
                                           int *displs_send_c, int *counts_recv_c,
                                           int *displs_recv_c, int *idx_recv_c, const char *varname_c)
{
    int ierr = MPI_SUCCESS;
    string varname = varname_c ? string(varname_c) : "var";

    // 将 Fortran 的 MPI_Fint 转为 C 的 MPI_Comm
    MPI_Comm all_comm = MPI_Comm_f2c((MPI_Fint)all_comm_f);

    int world_rank = -1;
    int world_size = -1;
    MPI_Comm_rank(all_comm, &world_rank);
    MPI_Comm_size(all_comm, &world_size);

    // 让每个相同Color进程拿到自己的rank_local
    MPI_Comm color_comm = MPI_COMM_NULL;
    ierr = MPI_Comm_split(all_comm, color_id, world_rank, &color_comm);
    if (ierr != MPI_SUCCESS)
        return;
    int rank_local = -1;
    MPI_Comm_rank(color_comm, &rank_local);
    MPI_Comm_free(&color_comm); // 先释放，后面重新创建

    // 根据rank_local重新创建新的组
    ierr = MPI_Comm_split(all_comm, rank_local, world_rank, &color_comm);
    // 获取新组的信息
    int nprocs_color = -1, comm_rank = -1;
    ierr = MPI_Comm_size(color_comm, &nprocs_color);
    ierr = MPI_Comm_rank(color_comm, &comm_rank);

    // 用一次 MPI_Alltoall 交换整张“行=src, 列=target”的矩阵
    const int *send_ptr = displs_send_c;
    vector<int> tmp_send;
    vector<int> displs_local(static_cast<size_t>(world_size), -1);
    if (send_ptr == nullptr)
    {
        tmp_send.assign(static_cast<size_t>(world_size), -1);
        send_ptr = tmp_send.data();
    }
    ierr = MPI_Alltoall(send_ptr, 1, MPI_INT,
                        displs_local.data(), 1, MPI_INT,
                        all_comm);
    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }

    // 首先收集每个进程的元素数
    vector<int> recvcounts, displs;
    vector<int> world_ranks; // 新增：收集所有进程的world_rank

    if (color_id == 0)
    {
        recvcounts.resize(nprocs_color, 0);
        displs.resize(nprocs_color, 0);
        world_ranks.resize(nprocs_color, 0);
    }

    int local_len_copy = local_len;
    ierr = MPI_Gather(&local_len_copy, 1, MPI_INT,
                      color_id == 0 ? recvcounts.data() : nullptr, 1, MPI_INT,
                      0, color_comm);

    // 收集所有进程的world_rank
    int world_rank_copy = world_rank;
    ierr = MPI_Gather(&world_rank_copy, 1, MPI_INT,
                      color_id == 0 ? world_ranks.data() : nullptr, 1, MPI_INT,
                      0, color_comm);

    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }

    int total_nums = 0;
    if (color_id == 0)
    {
        displs[0] = 0;
        for (int i = 1; i < nprocs_color; ++i)
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        total_nums = 0;
        for (int v : recvcounts)
            total_nums += v;
    }

    int send_elems = local_len;
    // 根进程接收缓冲
    vector<float> recv_bytes;
    if (color_id == 0)
    {
        recv_bytes.resize(static_cast<size_t>(total_nums));
    }
    // Gatherv: 单位是元素
    ierr = MPI_Gatherv(temp, send_elems, MPI_FLOAT,
                       color_id == 0 ? recv_bytes.data() : nullptr,
                       color_id == 0 ? recvcounts.data() : nullptr,
                       color_id == 0 ? displs.data() : nullptr,
                       MPI_FLOAT, 0, color_comm);
    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }

    if (color_id == 0)
    {
        string s_outdir = outdir ? outdir : "";
        string s_file = filename ? filename : "";
        // string path = s_outdir + s_file + "_" + varname + "_grp" + to_string(rank_local) + ".bin";
        string path = s_outdir + varname + ".dir/" +
                      s_file + "_" + varname + "_grp" + to_string(rank_local) + ".bin";

        ios::openmode mode = ios::binary | ios::out;
        mode |= (open_mode == 1 ? ios::trunc : ios::app);

        ofstream ofs(path, mode);
        if (!ofs)
        {
            MPI_Comm_free(&color_comm);
            return;
        }

        // 写入头信息：my_count, my_start, group_procs_num, total_nums
        int32_t mycount = static_cast<int32_t>(my_count);
        int32_t mystart = static_cast<int32_t>(my_start);
        int32_t world_procs = static_cast<int32_t>(world_size);
        int32_t group_procs = static_cast<int32_t>(nprocs_color);
        int32_t totalNums = static_cast<int32_t>(total_nums);

        ofs.write(reinterpret_cast<const char *>(&mycount), sizeof(mycount));
        ofs.write(reinterpret_cast<const char *>(&mystart), sizeof(mystart));
        ofs.write(reinterpret_cast<const char *>(&world_procs), sizeof(world_procs));
        ofs.write(reinterpret_cast<const char *>(&group_procs), sizeof(group_procs));
        ofs.write(reinterpret_cast<const char *>(&totalNums), sizeof(totalNums));

        // 写入列索引数组 counts_recv_c, displs_recv_c, idx_recv_c, displs_local
        // 大小分别是 world_size, world_size, my_count
        if (world_size > 0)
        {
            ofs.write(reinterpret_cast<const char *>(counts_recv_c),
                      static_cast<streamsize>(world_size * sizeof(int)));
            ofs.write(reinterpret_cast<const char *>(displs_recv_c),
                      static_cast<streamsize>(world_size * sizeof(int)));
            ofs.write(reinterpret_cast<const char *>(displs_local.data()),
                      static_cast<streamsize>(world_size * sizeof(int)));
        }
        if (my_count > 0)
        {
            ofs.write(reinterpret_cast<const char *>(idx_recv_c),
                      static_cast<streamsize>(my_count * sizeof(int)));
        }
        // 写入数据来源进程的 world_rank 和对应的列数
        if (group_procs > 0)
        {
            // world_ranks 和 recvcounts 已由 MPI_Gather 在 root 收集
            for (int i = 0; i < group_procs; ++i)
            {
                int32_t wr = static_cast<int32_t>(world_ranks[i]);
                int32_t lc = static_cast<int32_t>(recvcounts[i]);
                ofs.write(reinterpret_cast<const char *>(&wr), sizeof(wr));
                ofs.write(reinterpret_cast<const char *>(&lc), sizeof(lc));
            }
        }
        // 写数据（列主序的线性存储）
        if (!recv_bytes.empty())
            ofs.write(reinterpret_cast<const char *>(recv_bytes.data()),
                      static_cast<streamsize>(recv_bytes.size() * sizeof(float)));
        ofs.close();
    }
    MPI_Comm_free(&color_comm);
}

extern "C" void wrap_output_2d_gather_sp_c(const float *temp, int len_dim1, int local_cols,
                                           int color_id, int all_comm_f,
                                           const char *outdir, const char *filename,
                                           int open_mode,
                                           int my_count, int my_start,
                                           int *displs_send_c, int *counts_recv_c,
                                           int *displs_recv_c, int *idx_recv_c, const char *varname_c)
{
    int ierr = MPI_SUCCESS;
    string varname = varname_c ? string(varname_c) : "var";

    // 将 Fortran 的 MPI_Fint 转为 C 的 MPI_Comm
    MPI_Comm all_comm = MPI_Comm_f2c((MPI_Fint)all_comm_f);

    int world_rank = -1;
    int world_size = -1;
    MPI_Comm_rank(all_comm, &world_rank);
    MPI_Comm_size(all_comm, &world_size);

    // 让每个相同Color进程拿到自己的rank_local
    MPI_Comm color_comm = MPI_COMM_NULL;
    ierr = MPI_Comm_split(all_comm, color_id, world_rank, &color_comm);
    if (ierr != MPI_SUCCESS)
        return;
    int rank_local = -1;
    MPI_Comm_rank(color_comm, &rank_local);
    MPI_Comm_free(&color_comm); // 先释放，后面重新创建

    // 根据rank_local重新创建新的组
    ierr = MPI_Comm_split(all_comm, rank_local, world_rank, &color_comm);
    // 获取新组的信息
    int nprocs_color = -1;
    ierr = MPI_Comm_size(color_comm, &nprocs_color);

    // printf("world rank: %d, color_id: %d, displs_send[0]: %d\n",
    //        world_rank, color_id,
    //        displs_send_c ? displs_send_c[0] : -1);

    // 用一次 MPI_Alltoall 交换整张“行=src, 列=target”的矩阵
    const int *send_ptr = displs_send_c;
    vector<int> tmp_send;
    vector<int> displs_local(static_cast<size_t>(world_size), -1);
    if (send_ptr == nullptr)
    {
        tmp_send.assign(static_cast<size_t>(world_size), -1);
        send_ptr = tmp_send.data();
    }
    ierr = MPI_Alltoall(send_ptr, 1, MPI_INT,
                        displs_local.data(), 1, MPI_INT,
                        all_comm);
    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }

    // 每进程元素个数（float）
    const int send_elems = len_dim1 * local_cols;

    // 首先收集每个进程的列数
    vector<int> recvcounts_cols, displs_cols;
    vector<int> world_ranks; // 新增：收集所有进程的world_rank

    if (color_id == 0)
    {
        recvcounts_cols.resize(nprocs_color, 0);
        displs_cols.resize(nprocs_color, 0);
        world_ranks.resize(nprocs_color, 0);
        // printf("[DEBUG] Color ID %d nprocs_color: %d\n", color_id, nprocs_color);
    }

    int local_cols_copy = local_cols;
    ierr = MPI_Gather(&local_cols_copy, 1, MPI_INT,
                      color_id == 0 ? recvcounts_cols.data() : nullptr, 1, MPI_INT,
                      0, color_comm);

    // 收集所有进程的world_rank
    int world_rank_copy = world_rank;
    ierr = MPI_Gather(&world_rank_copy, 1, MPI_INT,
                      color_id == 0 ? world_ranks.data() : nullptr, 1, MPI_INT,
                      0, color_comm);

    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }

    int total_cols = 0;
    vector<int> recvcounts_elems, displs_elems;
    if (color_id == 0)
    {
        displs_cols[0] = 0;
        for (int i = 1; i < nprocs_color; ++i)
            displs_cols[i] = displs_cols[i - 1] + recvcounts_cols[i - 1];
        total_cols = 0;
        for (int v : recvcounts_cols)
            total_cols += v;

        recvcounts_elems.resize(nprocs_color, 0);
        displs_elems.resize(nprocs_color, 0);
        for (int i = 0; i < nprocs_color; ++i)
        {
            recvcounts_elems[i] = recvcounts_cols[i] * len_dim1;
            displs_elems[i] = displs_cols[i] * len_dim1;
        }
    }

    // 根进程接收缓冲
    vector<float> recvbuf;
    if (color_id == 0)
    {
        recvbuf.resize(static_cast<size_t>(len_dim1) * static_cast<size_t>(total_cols));
        // printf("[DEBUG] Color ID %d total_cols: %d, recvbuf size: %zu\n",
        //    color_id, total_cols, recvbuf.size());
    }

    // Gatherv: 单位是 float 元素
    ierr = MPI_Gatherv(temp, send_elems, MPI_FLOAT,
                       color_id == 0 ? recvbuf.data() : nullptr,
                       color_id == 0 ? recvcounts_elems.data() : nullptr,
                       color_id == 0 ? displs_elems.data() : nullptr,
                       MPI_FLOAT, 0, color_comm);

    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }

    // 根进程写文件：<outdir><filename>_grp<color_id>.bin
    if (color_id == 0)
    {
        string s_outdir = outdir ? outdir : "";
        string s_file = filename ? filename : "";
        // string path = s_outdir + s_file + "_" + varname + "_grp" + to_string(rank_local) + ".bin";
        string path = s_outdir + varname + ".dir/" +
                      s_file + "_" + varname + "_grp" + to_string(rank_local) + ".bin";
        // printf("[INFO] Writing binary file: %s\n", path.c_str());

        ios::openmode mode = ios::binary | ios::out;
        mode |= (open_mode == 1 ? ios::trunc : ios::app);

        ofstream ofs(path, mode);
        if (!ofs)
        {
            // 无法打开就直接返回（避免在并行中抛异常）
            MPI_Comm_free(&color_comm);
            return;
        }

        // 写入头信息：my_count, my_start, group_procs_num, total_cols, len_dim1
        int32_t mycount = static_cast<int32_t>(my_count);
        int32_t mystart = static_cast<int32_t>(my_start);
        int32_t world_procs_num = static_cast<int32_t>(world_size);
        int32_t group_procs_num = static_cast<int32_t>(nprocs_color);
        int32_t totalCols = static_cast<int32_t>(total_cols);
        int32_t dim1_lens = static_cast<int32_t>(len_dim1);
        ofs.write(reinterpret_cast<const char *>(&mycount), sizeof(mycount));
        ofs.write(reinterpret_cast<const char *>(&mystart), sizeof(mystart));
        ofs.write(reinterpret_cast<const char *>(&world_procs_num), sizeof(world_procs_num));
        ofs.write(reinterpret_cast<const char *>(&group_procs_num), sizeof(group_procs_num));
        ofs.write(reinterpret_cast<const char *>(&totalCols), sizeof(totalCols));
        ofs.write(reinterpret_cast<const char *>(&dim1_lens), sizeof(dim1_lens));

        // 写入列索引数组 counts_recv_c, displs_recv_c, idx_recv_c, displs_local
        // 大小分别是 world_size, world_size, my_count
        if (world_size > 0)
        {
            ofs.write(reinterpret_cast<const char *>(counts_recv_c),
                      static_cast<streamsize>(world_size * sizeof(int)));
            ofs.write(reinterpret_cast<const char *>(displs_recv_c),
                      static_cast<streamsize>(world_size * sizeof(int)));
            ofs.write(reinterpret_cast<const char *>(displs_local.data()),
                      static_cast<streamsize>(world_size * sizeof(int)));
        }
        if (my_count > 0)
        {
            ofs.write(reinterpret_cast<const char *>(idx_recv_c),
                      static_cast<streamsize>(my_count * sizeof(int)));
        }
        // 写入数据来源进程的 world_rank 和对应的列数
        if (group_procs_num > 0)
        {
            // world_ranks 和 recvcounts_cols 已由 MPI_Gather 在 root 收集
            for (int i = 0; i < group_procs_num; ++i)
            {
                int32_t wr = static_cast<int32_t>(world_ranks[i]);
                int32_t lc = static_cast<int32_t>(recvcounts_cols[i]);
                ofs.write(reinterpret_cast<const char *>(&wr), sizeof(wr));
                ofs.write(reinterpret_cast<const char *>(&lc), sizeof(lc));
            }
        }
        // 写数据（列主序的线性存储）
        if (!recvbuf.empty())
        {
            ofs.write(reinterpret_cast<const char *>(recvbuf.data()),
                      static_cast<streamsize>(recvbuf.size() * sizeof(float)));
        }
        ofs.close();
    }

    MPI_Comm_free(&color_comm);
}

extern "C" void wrap_output_3d_gather_sp_c(const float *temp, int len_dim1, int len_dim3,
                                           int local_cols, int color_id, int all_comm_f,
                                           const char *outdir, const char *filename,
                                           int open_mode,
                                           int my_count, int my_start,
                                           int *displs_send_c, int *counts_recv_c,
                                           int *displs_recv_c, int *idx_recv_c, const char *varname_c)
{
    int ierr = MPI_SUCCESS;
    string varname = varname_c ? string(varname_c) : "var";

    // 将 Fortran 的 MPI_Fint 转为 C 的 MPI_Comm
    MPI_Comm all_comm = MPI_Comm_f2c((MPI_Fint)all_comm_f);

    int world_rank = -1;
    int world_size = -1;
    MPI_Comm_rank(all_comm, &world_rank);
    MPI_Comm_size(all_comm, &world_size);

    // 让每个相同Color进程拿到自己的rank_local
    MPI_Comm color_comm = MPI_COMM_NULL;
    ierr = MPI_Comm_split(all_comm, color_id, world_rank, &color_comm);
    if (ierr != MPI_SUCCESS)
        return;
    int rank_local = -1;
    MPI_Comm_rank(color_comm, &rank_local);
    MPI_Comm_free(&color_comm); // 先释放，后面重新创建

    // 根据rank_local重新创建新的组
    ierr = MPI_Comm_split(all_comm, rank_local, world_rank, &color_comm);
    // 获取新组的信息
    int nprocs_color = -1;
    ierr = MPI_Comm_size(color_comm, &nprocs_color);

    // printf("world rank: %d, color_id: %d, displs_send[0]: %d\n",
    //        world_rank, color_id,
    //        displs_send_c ? displs_send_c[0] : -1);

    // 用一次 MPI_Alltoall 交换整张“行=src, 列=target”的矩阵
    const int *send_ptr = displs_send_c;
    vector<int> tmp_send;
    vector<int> displs_local(static_cast<size_t>(world_size), -1);
    if (send_ptr == nullptr)
    {
        tmp_send.assign(static_cast<size_t>(world_size), -1);
        send_ptr = tmp_send.data();
    }
    ierr = MPI_Alltoall(send_ptr, 1, MPI_INT,
                        displs_local.data(), 1, MPI_INT,
                        all_comm);
    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }
    
    // 每列元素个数（float）
    const int elem_per_col = len_dim1 * len_dim3;
    // 每进程元素个数（float）
    const int send_elems = local_cols * elem_per_col;

    // 首先收集每个进程的列数
    vector<int> recvcounts_cols, displs_cols;
    vector<int> world_ranks; // 新增：收集所有进程的world_rank

    if (color_id == 0)
    {
        recvcounts_cols.resize(nprocs_color, 0);
        displs_cols.resize(nprocs_color, 0);
        world_ranks.resize(nprocs_color, 0);
        // printf("[DEBUG] Color ID %d nprocs_color: %d\n", color_id, nprocs_color);
    }

    int local_cols_copy = local_cols;
    ierr = MPI_Gather(&local_cols_copy, 1, MPI_INT,
                      color_id == 0 ? recvcounts_cols.data() : nullptr, 1, MPI_INT,
                      0, color_comm);

    // 收集所有进程的world_rank
    int world_rank_copy = world_rank;
    ierr = MPI_Gather(&world_rank_copy, 1, MPI_INT,
                      color_id == 0 ? world_ranks.data() : nullptr, 1, MPI_INT,
                      0, color_comm);

    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }

    int total_cols = 0;
    vector<int> recvcounts_elems, displs_elems;
    if (color_id == 0)
    {
        displs_cols[0] = 0;
        for (int i = 1; i < nprocs_color; ++i)
            displs_cols[i] = displs_cols[i - 1] + recvcounts_cols[i - 1];
        total_cols = 0;
        for (int v : recvcounts_cols)
            total_cols += v;

        recvcounts_elems.resize(nprocs_color, 0);
        displs_elems.resize(nprocs_color, 0);
        for (int i = 0; i < nprocs_color; ++i)
        {
            recvcounts_elems[i] = recvcounts_cols[i] * elem_per_col;
            displs_elems[i] = displs_cols[i] * elem_per_col;
        }
    }

    // 根进程接收缓冲
    vector<float> recvbuf;
    if (color_id == 0)
    {
        recvbuf.resize(static_cast<size_t>(elem_per_col) * static_cast<size_t>(total_cols));
        // printf("[DEBUG] Color ID %d total_cols: %d, recvbuf size: %zu\n",
        //    color_id, total_cols, recvbuf.size());
    }

    // Gatherv: 单位是 float 元素
    ierr = MPI_Gatherv(temp, send_elems, MPI_FLOAT,
                       color_id == 0 ? recvbuf.data() : nullptr,
                       color_id == 0 ? recvcounts_elems.data() : nullptr,
                       color_id == 0 ? displs_elems.data() : nullptr,
                       MPI_FLOAT, 0, color_comm);

    if (ierr != MPI_SUCCESS)
    {
        MPI_Comm_free(&color_comm);
        return;
    }

    // 根进程写文件：<outdir><filename>_grp<color_id>.bin
    if (color_id == 0)
    {
        string s_outdir = outdir ? outdir : "";
        string s_file = filename ? filename : "";
        // string path = s_outdir + s_file + "_" + varname + "_grp" + to_string(rank_local) + ".bin";
        string path = s_outdir + varname + ".dir/" +
                      s_file + "_" + varname + "_grp" + to_string(rank_local) + ".bin";
        // printf("[INFO] Writing binary file: %s\n", path.c_str());

        ios::openmode mode = ios::binary | ios::out;
        mode |= (open_mode == 1 ? ios::trunc : ios::app);

        ofstream ofs(path, mode);
        if (!ofs)
        {
            // 无法打开就直接返回（避免在并行中抛异常）
            MPI_Comm_free(&color_comm);
            return;
        }

        // 写入头信息：my_count, my_start, group_procs_num, total_cols, len_dim1
        int32_t mycount = static_cast<int32_t>(my_count);
        int32_t mystart = static_cast<int32_t>(my_start);
        int32_t world_procs_num = static_cast<int32_t>(world_size);
        int32_t group_procs_num = static_cast<int32_t>(nprocs_color);
        int32_t totalCols = static_cast<int32_t>(total_cols);
        int32_t dim1_lens = static_cast<int32_t>(len_dim1);
        int32_t dim3_lens = static_cast<int32_t>(len_dim3);

        ofs.write(reinterpret_cast<const char *>(&mycount), sizeof(mycount));
        ofs.write(reinterpret_cast<const char *>(&mystart), sizeof(mystart));
        ofs.write(reinterpret_cast<const char *>(&world_procs_num), sizeof(world_procs_num));
        ofs.write(reinterpret_cast<const char *>(&group_procs_num), sizeof(group_procs_num));
        ofs.write(reinterpret_cast<const char *>(&totalCols), sizeof(totalCols));
        ofs.write(reinterpret_cast<const char *>(&dim1_lens), sizeof(dim1_lens));
        ofs.write(reinterpret_cast<const char *>(&dim3_lens), sizeof(dim3_lens));

        // 写入列索引数组 counts_recv_c, displs_recv_c, idx_recv_c, displs_local
        // 大小分别是 world_size, world_size, my_count
        if (world_size > 0)
        {
            ofs.write(reinterpret_cast<const char *>(counts_recv_c),
                      static_cast<streamsize>(world_size * sizeof(int)));
            ofs.write(reinterpret_cast<const char *>(displs_recv_c),
                      static_cast<streamsize>(world_size * sizeof(int)));
            ofs.write(reinterpret_cast<const char *>(displs_local.data()),
                      static_cast<streamsize>(world_size * sizeof(int)));
        }
        if (my_count > 0)
        {
            ofs.write(reinterpret_cast<const char *>(idx_recv_c),
                      static_cast<streamsize>(my_count * sizeof(int)));
        }
        // 写入数据来源进程的 world_rank 和对应的列数
        if (group_procs_num > 0)
        {
            // world_ranks 和 recvcounts_cols 已由 MPI_Gather 在 root 收集
            for (int i = 0; i < group_procs_num; ++i)
            {
                int32_t wr = static_cast<int32_t>(world_ranks[i]);
                int32_t lc = static_cast<int32_t>(recvcounts_cols[i]);
                ofs.write(reinterpret_cast<const char *>(&wr), sizeof(wr));
                ofs.write(reinterpret_cast<const char *>(&lc), sizeof(lc));
            }
        }
        // 写数据（列主序：每列 elem_per_col 个 float）
        if (!recvbuf.empty())
        {
            ofs.write(reinterpret_cast<const char *>(recvbuf.data()),
                      static_cast<streamsize>(recvbuf.size() * sizeof(float)));
        }
        ofs.close();
    }

    MPI_Comm_free(&color_comm);
}

extern "C" void bin_to_nc_group_sp(const char *outdir, int all_comm_f, int color, int nvars, const char *varnames_c)
{
    int ierr = MPI_SUCCESS;

    // 将 Fortran 的 MPI_Fint 转为 C 的 MPI_Comm
    MPI_Comm all_comm = MPI_Comm_f2c((MPI_Fint)all_comm_f);

    int world_size = -1;
    int world_rank = -1;
    MPI_Comm_size(all_comm, &world_size);
    MPI_Comm_rank(all_comm, &world_rank);

    MPI_Comm color_comm = MPI_COMM_NULL;
    ierr = MPI_Comm_split(all_comm, color, world_rank, &color_comm);

    if (ierr != MPI_SUCCESS)
    {
        printf("[ERROR] MPI_Comm_split failed: %d\n", ierr);
        return;
    }

    int rank_local = -1;
    int nprocs_color = -1;
    MPI_Comm_rank(color_comm, &rank_local);
    MPI_Comm_size(color_comm, &nprocs_color);

    if (color == 0)
    {
        for (int vars_i = 0; vars_i < nvars; vars_i++)
        {
            string varname;
            if (varnames_c != nullptr)
            {
                size_t varname_len = strlen(varnames_c);
                size_t offset = static_cast<size_t>(vars_i) * (varname_len + 1);
                varname = string(varnames_c + offset);
            }
            string s_outdir = outdir ? outdir : "";
            string bin_path = s_outdir + "_" + varname + "_grp" + to_string(rank_local) + ".bin";
            ifstream ifs(bin_path, ios::binary);
            if (!ifs)
            {
                printf("[WARN] Cannot open bin file: %s (skip)\n", bin_path.c_str());
                MPI_Comm_free(&color_comm);
                return;
            }
            printf("[INFO] Processing variable '%s' from file: %s\n", varname.c_str(), bin_path.c_str());
            // 读取头信息（按写入时的顺序）
            int32_t mycount = 0, mystart = 0, world_procs_num = 0, group_procs_num = 0, totalCols = 0, dim1_lens = 0;
            ifs.read(reinterpret_cast<char *>(&mycount), sizeof(mycount));
            ifs.read(reinterpret_cast<char *>(&mystart), sizeof(mystart));
            ifs.read(reinterpret_cast<char *>(&world_procs_num), sizeof(world_procs_num));
            ifs.read(reinterpret_cast<char *>(&group_procs_num), sizeof(group_procs_num));
            ifs.read(reinterpret_cast<char *>(&totalCols), sizeof(totalCols));
            ifs.read(reinterpret_cast<char *>(&dim1_lens), sizeof(dim1_lens));
            if (ifs.fail())
            {
                printf("[ERROR] Failed to read header from %s\n", bin_path.c_str());
                ifs.close();
                MPI_Comm_free(&color_comm);
                return;
            }

            // 读取 counts_recv_c 、 displs_recv_c 、displs_local （长度为 world_size）
            vector<int> counts_recv(world_procs_num), displs_recv(world_procs_num), displs_local(world_procs_num);
            if (world_procs_num > 0)
            {
                ifs.read(reinterpret_cast<char *>(counts_recv.data()), static_cast<streamsize>(world_procs_num * sizeof(int)));
                ifs.read(reinterpret_cast<char *>(displs_recv.data()), static_cast<streamsize>(world_procs_num * sizeof(int)));
                ifs.read(reinterpret_cast<char *>(displs_local.data()), static_cast<streamsize>(world_procs_num * sizeof(int)));
                if (ifs.fail())
                {
                    printf("[ERROR] Failed to read counts/displs from %s\n", bin_path.c_str());
                    ifs.close();
                }
            }

            // 读取 idx_recv_c（长度为 mycount）
            vector<int> idx_recv;
            if (!ifs.fail() && mycount > 0)
            {
                idx_recv.resize(static_cast<size_t>(mycount));
                ifs.read(reinterpret_cast<char *>(idx_recv.data()), static_cast<streamsize>(mycount * sizeof(int)));
                if (ifs.fail())
                {
                    printf("[ERROR] Failed to read idx_recv from %s\n", bin_path.c_str());
                    ifs.close();
                }
            }

            // 读取来源进程信息（group_procs_num 个二元组：world_rank, local_cols）
            vector<int32_t> src_world_ranks(static_cast<size_t>(group_procs_num));
            vector<int32_t> src_recv_count(static_cast<size_t>(group_procs_num));
            if (!ifs.fail())
            {
                for (int i = 0; i < group_procs_num; ++i)
                {
                    ifs.read(reinterpret_cast<char *>(&src_world_ranks[i]), sizeof(int32_t));
                    ifs.read(reinterpret_cast<char *>(&src_recv_count[i]), sizeof(int32_t));
                    if (ifs.fail())
                    {
                        printf("[ERROR] Failed to read source info from %s\n", bin_path.c_str());
                        break;
                    }
                }
            }

            // 元素数量可由 src_recv_count（每个 bin 中记录的列数）和 dim1_lens 计算得到
            vector<float *> data_ptrs(static_cast<size_t>(world_procs_num), nullptr);
            // 读取数据区：dim1_lens * totalCols 个 float（列主序）
            if (!ifs.fail())
            {
                bool read_error = false;
                for (int i = 0; i < group_procs_num; ++i)
                {
                    int wr = static_cast<int>(src_world_ranks[i]);
                    int cols = static_cast<int>(src_recv_count[i]);
                    size_t nelems_i = static_cast<size_t>(dim1_lens) * static_cast<size_t>(cols);
                    if (nelems_i > 0)
                    {
                        float *buf = new (nothrow) float[nelems_i];
                        if (buf == nullptr)
                        {
                            printf("[ERROR] Out of memory allocating %zu floats for world_rank %d\n", nelems_i, wr);
                            read_error = true;
                            break;
                        }
                        ifs.read(reinterpret_cast<char *>(buf), static_cast<streamsize>(nelems_i * sizeof(float)));
                        if (ifs.fail())
                        {
                            printf("[ERROR] Failed to read data block for world_rank %d from %s\n", wr, bin_path.c_str());
                            delete[] buf;
                            read_error = true;
                            break;
                        }
                        // 存储指针（后续若需元素数量，可用 src_recv_count[对应index] * dim1_lens 计算）
                        data_ptrs[wr] = buf;
                    }
                    // cols == 0 时保持 nullptr
                }
                ifs.close();
            }

            // if (rank_local == 1)
            // {
            //     // for(int i = 0; i < group_procs_num; i++)
            //     // {
            //     //     printf("Debug: src_world_ranks[%d] = %d, src_recv_count[%d] = %d\n",
            //     //            i, src_world_ranks[i],
            //     //            i, src_recv_count[i]);
            //     // }
            //     printf("Debug: counts_recv for world procs: %d\n", src_recv_count[1]);
            //     float *ptr = data_ptrs[16];
            //     for (int i = 0; i < src_recv_count[1]; i++)
            //     {
            //         for (int j = 0; j < dim1_lens; j++)
            //         {
            //             printf("data_ptrs[16][%d][%d] = %f\n", i, j, ptr[i * dim1_lens + j]);
            //         }
            //     }
            // }

            // 交换 world_rank 拥有数据的情况
            vector<int> owner_has(static_cast<size_t>(world_procs_num), 0);
            vector<int> all_owner_has(static_cast<size_t>(nprocs_color) * static_cast<size_t>(world_procs_num), 0);
            vector<int> owner_of_world(static_cast<size_t>(world_procs_num), -1);
            for (int i = 0; i < group_procs_num; ++i)
            {
                int wr = static_cast<int>(src_world_ranks[i]);
                if (wr >= 0 && wr < world_procs_num)
                    owner_has[static_cast<size_t>(wr)] = 1; // 若同一 wr 出现多次，仍设为1
            }
            MPI_Allgather(owner_has.data(), world_procs_num, MPI_INT,
                          all_owner_has.data(), world_procs_num, MPI_INT,
                          color_comm);
            // 构建 world_rank -> color0_local_rank 映射（若多个 owner 拥有同一 world_rank，则警告）
            for (int owner = 0; owner < nprocs_color; ++owner)
            {
                for (int wr = 0; wr < world_procs_num; ++wr)
                {
                    int has = all_owner_has[owner * world_procs_num + wr];
                    if (has)
                    {
                        if (owner_of_world[wr] == -1)
                            owner_of_world[wr] = owner;
                        else if (owner_of_world[wr] != owner)
                            printf("[WARN] world_rank %d found in multiple color0 owners: first=%d other=%d\n",
                                   wr, owner_of_world[wr], owner);
                    }
                }
            }

            vector<float> recv_buffer;
            if (mycount > 0)
            {
                recv_buffer.resize(static_cast<size_t>(dim1_lens) * static_cast<size_t>(mycount));
            }

            // 发起接收：遍历所有 world_rank，若 counts_recv[wr] > 0，则从对应 owner 接收子块
            vector<MPI_Request> recv_reqs;
            for (int wr = 0; wr < world_procs_num; ++wr)
            {
                int need_cols = counts_recv[wr];
                if (need_cols <= 0)
                    continue;

                int owner = owner_of_world[wr];
                if (owner == -1)
                {
                    printf("[ERROR] No owner found for world_rank %d\n", wr);
                    continue;
                }

                // 计算接收位置：从 displs_recv[wr] 列开始，接收 need_cols 列
                size_t recv_col_start = static_cast<size_t>(displs_recv[wr]);
                size_t recv_elem_start = recv_col_start * static_cast<size_t>(dim1_lens);
                size_t recv_nelems = static_cast<size_t>(need_cols) * static_cast<size_t>(dim1_lens);

                if (recv_elem_start + recv_nelems > recv_buffer.size())
                {
                    printf("[ERROR] Receive buffer overflow for wr %d\n", wr);
                    continue;
                }

                MPI_Request req;
                MPI_Irecv(recv_buffer.data() + recv_elem_start, static_cast<int>(recv_nelems), MPI_FLOAT,
                          owner, /*tag=*/wr, color_comm, &req);
                recv_reqs.push_back(req);
            }

            // 从 all_owner_has 推断 target 是否需要 wr（但 all_owner_has 是拥有情况，不是需要情况）
            // 实际需要从 all_counts_recv[target * world_size + wr] 获取，但我们没有 all_counts_recv
            // 先收集所有 counts_recv
            // 发起发送：遍历本进程持有的 wr，向需要该 wr 的其他进程发送子块

            vector<int> all_displs_local(nprocs_color * world_procs_num, 0);
            MPI_Allgather(displs_local.data(), world_procs_num, MPI_INT,
                          all_displs_local.data(), world_procs_num, MPI_INT,
                          color_comm);

            vector<int> all_counts_recv(static_cast<size_t>(nprocs_color) * static_cast<size_t>(world_procs_num), 0);
            MPI_Allgather(counts_recv.data(), world_procs_num, MPI_INT,
                          all_counts_recv.data(), world_procs_num, MPI_INT,
                          color_comm);

            vector<MPI_Request> send_reqs;
            for (int wr = 0; wr < world_procs_num; ++wr)
            {
                float *data_block = data_ptrs[wr];
                if (data_block == nullptr)
                    continue;

                // 找到该 wr 的所有需要者（通过反向查询 all_owner_has 或遍历所有进程）
                // 这里简化：遍历所有 color0 进程，检查其 counts_recv[wr] > 0
                for (int target = 0; target < nprocs_color; ++target)
                {
                    int target_need = all_counts_recv[target * world_procs_num + wr];
                    int target_start = all_displs_local[target * world_procs_num + wr];
                    // 计算发送子块：从 displs_local[wr] 列开始，发送 target_need 列
                    if (target_need <= 0)
                        continue;

                    size_t send_elem_start = static_cast<size_t>(target_start) * static_cast<size_t>(dim1_lens);
                    size_t send_nelems = static_cast<size_t>(target_need) * static_cast<size_t>(dim1_lens);

                    MPI_Request req;
                    MPI_Isend(data_block + send_elem_start, static_cast<int>(send_nelems), MPI_FLOAT,
                              target, /*tag=*/wr, color_comm, &req);
                    send_reqs.push_back(req);
                }
            }

            // 等待所有通信完成
            if (!recv_reqs.empty())
                MPI_Waitall(static_cast<int>(recv_reqs.size()), recv_reqs.data(), MPI_STATUSES_IGNORE);
            if (!send_reqs.empty())
                MPI_Waitall(static_cast<int>(send_reqs.size()), send_reqs.data(), MPI_STATUSES_IGNORE);

            // 此时 recv_buffer 已填充完整，可根据 idx_recv 排序
            vector<float> put_buffer;
            put_buffer.resize((size_t)mycount * (size_t)dim1_lens);
            // reorder
            for (int i = 0; i < mycount; ++i)
            {
                int lidx = idx_recv[i] - 1; // 这一列的全局列号
                if (lidx < 0 || lidx >= mycount)
                {
                    fprintf(stderr, "[ERROR] idx_recv[%d] = %d out of range [0, %d)\n", i, lidx, mycount);
                    // MPI_Abort(color_comm, 1); // 立即停止，避免后续更复杂的内存破坏
                }
                else
                {
                    // src in recv_buffer
                    const float *src = recv_buffer.data() + (size_t)i * (size_t)dim1_lens;
                    // dest in put_buffer
                    float *dest = put_buffer.data() + (size_t)lidx * (size_t)dim1_lens;
                    memcpy(dest, src, (size_t)dim1_lens * sizeof(float));
                }
            }

            // 写 NetCDF 文件
            int ncid;
            int ierr;
            MPI_Offset start[2] = {static_cast<MPI_Offset>(mystart - 1), 0};
            MPI_Offset count[2] = {static_cast<MPI_Offset>(mycount), static_cast<MPI_Offset>(dim1_lens)};

            string nc_file = s_outdir + "_" + varname + ".nc";

            ierr = ncmpi_open(color_comm, nc_file.c_str(), NC_WRITE, MPI_INFO_NULL, &ncid);
            if (ierr != NC_NOERR)
            {
                printf("nc_open_par error: %s\n", ncmpi_strerror(ierr));
            }
            int nc_varid;
            ierr = ncmpi_inq_varid(ncid, varname.c_str(), &nc_varid);

            ierr = ncmpi_put_vara_float_all(ncid, nc_varid, start, count, put_buffer.data());
            if (ierr != NC_NOERR)
            {
                printf("PnetCDF error: %s\n", ncmpi_strerror(ierr));
            }

            ncmpi_close(ncid);
        }
    }
    MPI_Comm_free(&color_comm);
}