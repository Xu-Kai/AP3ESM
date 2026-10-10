#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <pnetcdf.h>
#include <mpi.h>
#include <algorithm>
#include <regex>

// #define FILE_NUM 16 // 每个进程最多处理20个二进制文件
int FILE_NUM = 16;

using namespace std;

struct BinFileData
{
    int bin_index;
    int32_t mycount, mystart, world_procs_num, group_procs_num, totalCols, dim1_lens;
    vector<int> counts_recv, displs_recv, displs_local;
    vector<int> idx_recv;
    vector<int32_t> src_world_ranks, src_recv_count;
    vector<float *> data_ptrs; // 按 world_rank 索引
};

void bin_to_nc_group_sp(string timestep, const char *outdir, int all_comm_f, int color, int nvars, const char *varnames_c)
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
                const char *p = varnames_c;
                for (int i = 0; i < vars_i; ++i)
                {
                    // 跳到下一个\0后面
                    p += strlen(p) + 1;
                }
                varname = string(p);
            }
            cout << "[INFO] Process " << world_rank << " handling variable '" << varname << "'\n";

            // 创建进程要读取的所有Bin文件的数据、信息
            vector<BinFileData> bin_files_data;
            int world_procs_num = 0;
            int group_procs_num = 0;

            for (int file_idx = 0; file_idx < FILE_NUM; file_idx++)
            {
                BinFileData current_file_data;
                int file_global_idx = rank_local * FILE_NUM + file_idx;
                current_file_data.bin_index = file_global_idx;

                string s_outdir = outdir ? outdir : "";
                // string bin_path = s_outdir + "_" + varname + "_grp" + to_string(file_global_idx) + ".bin";
                string bin_path = "./"  + varname + ".dir/" +
                                  s_outdir + "_" + varname + "_grp" +
                                  to_string(file_global_idx) + ".bin";

                // printf("[INFO] Processing bin file: %s\n", bin_path.c_str());

                ifstream ifs(bin_path, ios::binary);
                if (!ifs)
                {
                    // 文件是连续编号的，遇到不存在就停止
                    printf("[WARN] Cannot open bin file: %s (skip)\n", bin_path.c_str());
                    ifs.close();
                    break;
                }

                // printf("[INFO] Processing variable '%s' from file: %s\n", varname.c_str(), bin_path.c_str());
                // 读取头信息（按写入时的顺序）
                // int32_t mycount = 0, mystart = 0, world_procs_num = 0, group_procs_num = 0, totalCols = 0, dim1_lens = 0;
                ifs.read(reinterpret_cast<char *>(&current_file_data.mycount), sizeof(current_file_data.mycount));
                ifs.read(reinterpret_cast<char *>(&current_file_data.mystart), sizeof(current_file_data.mystart));
                ifs.read(reinterpret_cast<char *>(&current_file_data.world_procs_num), sizeof(current_file_data.world_procs_num));
                ifs.read(reinterpret_cast<char *>(&current_file_data.group_procs_num), sizeof(current_file_data.group_procs_num));
                ifs.read(reinterpret_cast<char *>(&current_file_data.totalCols), sizeof(current_file_data.totalCols));
                // ifs.read(reinterpret_cast<char *>(&current_file_data.dim1_lens), sizeof(current_file_data.dim1_lens));
                // 1D数据时，dim1_lens固定为1
                current_file_data.dim1_lens = 1;
                if (ifs.fail())
                {
                    printf("[ERROR] Failed to read header from %s\n", bin_path.c_str());
                    ifs.close();
                    MPI_Comm_free(&color_comm);
                    return;
                }
                // cout << "[INFO] Read header: BinIndex=" << file_global_idx
                //      << " mycount=" << current_file_data.mycount
                //      << " mystart=" << current_file_data.mystart
                //      << " world_procs_num=" << current_file_data.world_procs_num
                //      << " group_procs_num=" << current_file_data.group_procs_num
                //      << " totalCols=" << current_file_data.totalCols << "\n";
                // 读取 counts_recv_c 、 displs_recv_c 、displs_local （长度为 world_size）
                // vector<int> counts_recv(world_procs_num), displs_recv(world_procs_num), displs_local(world_procs_num);
                world_procs_num = current_file_data.world_procs_num;
                if (world_procs_num > 0)
                {
                    current_file_data.counts_recv.resize(world_procs_num);
                    current_file_data.displs_recv.resize(world_procs_num);
                    current_file_data.displs_local.resize(world_procs_num);

                    ifs.read(reinterpret_cast<char *>(current_file_data.counts_recv.data()), static_cast<streamsize>(world_procs_num * sizeof(int)));
                    ifs.read(reinterpret_cast<char *>(current_file_data.displs_recv.data()), static_cast<streamsize>(world_procs_num * sizeof(int)));
                    ifs.read(reinterpret_cast<char *>(current_file_data.displs_local.data()), static_cast<streamsize>(world_procs_num * sizeof(int)));
                    if (ifs.fail())
                    {
                        printf("[ERROR] Failed to read counts/displs from %s\n", bin_path.c_str());
                        ifs.close();
                    }
                }
                // cout << "[INFO] Read counts_recv, displs_recv, displs_local for " << world_procs_num << " procs.\n";
                // 读取 idx_recv_c（长度为 mycount）
                // vector<int> idx_recv;
                if (!ifs.fail() && current_file_data.mycount > 0)
                {
                    current_file_data.idx_recv.resize(static_cast<size_t>(current_file_data.mycount));
                    ifs.read(reinterpret_cast<char *>(current_file_data.idx_recv.data()),
                             static_cast<streamsize>(current_file_data.mycount * sizeof(int)));
                    if (ifs.fail())
                    {
                        printf("[ERROR] Failed to read idx_recv from %s\n", bin_path.c_str());
                        ifs.close();
                    }
                }

                // cout << "[INFO] Read idx_recv for " << current_file_data.mycount << " columns.\n";
                // 读取来源进程信息（group_procs_num 个二元组：world_rank, local_cols）
                group_procs_num = current_file_data.group_procs_num;
                // vector<int32_t> src_world_ranks(static_cast<size_t>(group_procs_num));
                // vector<int32_t> src_recv_count(static_cast<size_t>(group_procs_num));
                if (!ifs.fail())
                {
                    current_file_data.src_world_ranks.resize(static_cast<size_t>(group_procs_num));
                    current_file_data.src_recv_count.resize(static_cast<size_t>(group_procs_num));
                    for (int i = 0; i < group_procs_num; ++i)
                    {
                        ifs.read(reinterpret_cast<char *>(&current_file_data.src_world_ranks[i]), sizeof(int32_t));
                        ifs.read(reinterpret_cast<char *>(&current_file_data.src_recv_count[i]), sizeof(int32_t));
                        // cout << "[INFO] Read source info: binIndex: " << to_string(file_global_idx)
                        //      << " world_rank=" << current_file_data.src_world_ranks[i]
                        //      << " recv_count=" << current_file_data.src_recv_count[i] << "\n";
                        if (ifs.fail())
                        {
                            printf("[ERROR] Failed to read source info from %s\n", bin_path.c_str());
                            break;
                        }
                    }
                }
                // cout << "[INFO] Read source ranks and recv counts for " << group_procs_num << " group procs.\n";
                // 元素数量可由 src_recv_count（每个 bin 中记录的列数）和 dim1_lens 计算得到
                // vector<float *> data_ptrs(static_cast<size_t>(world_procs_num), nullptr);
                // 读取数据区：dim1_lens * totalCols 个 float（列主序）
                if (!ifs.fail())
                {
                    // 先分配空间，防止越界
                    current_file_data.data_ptrs.resize(static_cast<size_t>(world_procs_num), nullptr);

                    bool read_error = false;
                    for (int i = 0; i < group_procs_num; ++i)
                    {
                        int wr = static_cast<int>(current_file_data.src_world_ranks[i]);
                        int cols = static_cast<int>(current_file_data.src_recv_count[i]);
                        size_t nelems_i = static_cast<size_t>(current_file_data.dim1_lens) * static_cast<size_t>(cols);

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
                            current_file_data.data_ptrs[wr] = buf;
                        }
                    }
                    // cout << "[INFO] Read data blocks for " << group_procs_num << " source procs.\n";
                    bin_files_data.push_back(std::move(current_file_data));
                    ifs.close();
                }
            }

            // cout << "[INFO] Process " << world_rank << " finished reading bin files for variable '" << varname << "'\n";
            // 交换 world_rank 拥有数据的情况
            vector<int> owner_has(static_cast<size_t>(world_procs_num), 0);
            vector<int> all_owner_has(static_cast<size_t>(nprocs_color) * static_cast<size_t>(world_procs_num), 0);
            vector<int> owner_of_world(static_cast<size_t>(world_procs_num), -1);

            // 标记本进程拥有的 world_ranks
            for (const auto &file : bin_files_data)
            {
                for (size_t i = 0; i < file.src_world_ranks.size(); ++i)
                {
                    int wr = static_cast<int>(file.src_world_ranks[i]);
                    if (wr >= 0 && wr < world_procs_num)
                        owner_has[wr] = 1;
                }
            }

            // if (world_rank == 0)
            // {
            //     for (int wr = 0; wr < world_procs_num; ++wr)
            //     {
            //         cout << "[DEBUG] Process " << world_rank << " owner_has[" << wr << "] = " << owner_has[wr] << "\n";
            //     }
            // }

            // 所有进程拿到这个表
            MPI_Allgather(owner_has.data(), world_procs_num, MPI_INT,
                          all_owner_has.data(), world_procs_num, MPI_INT,
                          color_comm);

            // 构建 world_rank -> color0_local_rank 映射
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

            // if (world_rank == 0)
            // {
            //     for (int wr = 0; wr < world_procs_num; ++wr)
            //     {
            //         cout << "[DEBUG] Process " << world_rank << " owner_of_world[" << wr << "] = " << owner_of_world[wr] << "\n";
            //     }
            // }
            // cout << "[INFO] Process " << world_rank << " constructed owner_of_world mapping.\n";
            // vector<float> recv_buffer;
            vector<vector<float>> recv_buffers(bin_files_data.size());
            vector<MPI_Request> recv_reqs;
            for (int file_idx = 0; file_idx < bin_files_data.size(); ++file_idx)
            {
                BinFileData &file = bin_files_data[file_idx];
                if (file.mycount > 0)
                {
                    recv_buffers[file_idx].resize(
                        static_cast<size_t>(file.dim1_lens) *
                        static_cast<size_t>(file.mycount));
                }
                // 发起接收：遍历所有 world_rank，若 counts_recv[wr] > 0，则从对应 owner 接收子块
                for (int wr = 0; wr < file.world_procs_num; ++wr)
                {
                    int need_cols = file.counts_recv[wr];
                    if (need_cols <= 0)
                        continue;

                    int owner = owner_of_world[wr];
                    if (owner == -1)
                    {
                        printf("[ERROR] No owner found for world_rank %d\n", wr);
                        continue;
                    }

                    // 计算接收位置：从 displs_recv[wr] 列开始，接收 need_cols 列
                    size_t recv_col_start = static_cast<size_t>(file.displs_recv[wr]);
                    size_t recv_elem_start = recv_col_start * static_cast<size_t>(file.dim1_lens);
                    size_t recv_nelems = static_cast<size_t>(need_cols) * static_cast<size_t>(file.dim1_lens);

                    if (recv_elem_start + recv_nelems > recv_buffers[file_idx].size())
                    {
                        printf("[ERROR] Receive buffer overflow for wr %d\n", wr);
                        continue;
                    }

                    // 使用 world_rank 作为 tag，确保唯一性
                    int unique_tag = file_idx * file.world_procs_num + wr;

                    MPI_Request req;
                    MPI_Irecv(recv_buffers[file_idx].data() + recv_elem_start,
                              static_cast<int>(recv_nelems),
                              MPI_FLOAT,
                              owner,
                              unique_tag, // tags
                              color_comm,
                              &req);
                    recv_reqs.push_back(req);
                    // }
                }
            }

            // cout << "[INFO] Process " << world_rank << " posted all Irecv requests for variable '" << varname << "'\n";
            // 收集每个进程拥有的 bin 文件数量
            int local_bin_count = static_cast<int>(bin_files_data.size());
            vector<int> all_bin_counts(nprocs_color, 0);
            MPI_Allgather(&local_bin_count, 1, MPI_INT,
                          all_bin_counts.data(), 1, MPI_INT,
                          color_comm);

            vector<MPI_Request> send_reqs;
            vector<vector<int>> all_counts_recv_bin(FILE_NUM);
            vector<vector<int>> all_displs_local_bin(FILE_NUM);

            for (int i = 0; i < FILE_NUM; ++i)
            {
                if (i < bin_files_data.size())
                {
                    all_counts_recv_bin[i] = bin_files_data[i].counts_recv;
                    all_displs_local_bin[i] = bin_files_data[i].displs_local;
                }
                else
                {
                    // 补零，保证长度一致
                    all_counts_recv_bin[i] = vector<int>(world_procs_num, 0);
                    all_displs_local_bin[i] = vector<int>(world_procs_num, 0);
                }
            }
            // cout << "[INFO] Process " << world_rank << " prepared counts and displs for variable '" << varname << "'\n";

            // 展开为一维数组用于 MPI_Allgather
            vector<int> counts_recv_flat;
            vector<int> displs_local_flat;
            for (int i = 0; i < FILE_NUM; ++i)
            {
                counts_recv_flat.insert(counts_recv_flat.end(),
                                        all_counts_recv_bin[i].begin(),
                                        all_counts_recv_bin[i].end());
                displs_local_flat.insert(displs_local_flat.end(),
                                         all_displs_local_bin[i].begin(),
                                         all_displs_local_bin[i].end());
            }

            // 全局收集所有进程的所有 bin 文件信息
            vector<int> all_counts_recv_flat(nprocs_color * FILE_NUM * world_procs_num, 0);
            MPI_Allgather(counts_recv_flat.data(), static_cast<int>(counts_recv_flat.size()), MPI_INT,
                          all_counts_recv_flat.data(), static_cast<int>(counts_recv_flat.size()), MPI_INT,
                          color_comm);

            vector<int> all_displs_local_flat(nprocs_color * FILE_NUM * world_procs_num, 0);
            MPI_Allgather(displs_local_flat.data(), static_cast<int>(displs_local_flat.size()), MPI_INT,
                          all_displs_local_flat.data(), static_cast<int>(displs_local_flat.size()), MPI_INT,
                          color_comm);

            // cout << "[INFO] Process " << world_rank << " gathered all counts and displs for variable '" << varname << "'\n";
            // 当前拥有的bin文件
            for (int file_idx = 0; file_idx < bin_files_data.size(); ++file_idx)
            {
                BinFileData &file = bin_files_data[file_idx];
                // 遍历所有 world_rank，查看哪些进程需要数据
                for (int wr = 0; wr < file.world_procs_num; ++wr)
                {
                    float *data_block = file.data_ptrs[wr];
                    if (data_block == nullptr)
                        continue;
                    // 遍历所有目标进程
                    for (int target = 0; target < nprocs_color; ++target)
                    {
                        // 遍历所有目标进程的 bin 文件
                        for (int t_file_idx = 0; t_file_idx < all_bin_counts[target]; ++t_file_idx)
                        {
                            // 新的全局索引
                            int flat_idx = target * FILE_NUM * world_procs_num +
                                           t_file_idx * file.world_procs_num +
                                           wr;

                            int target_need = all_counts_recv_flat[flat_idx];
                            int target_start = all_displs_local_flat[flat_idx];

                            if (target_need <= 0)
                                continue;

                            size_t send_elem_start = static_cast<size_t>(target_start) *
                                                     static_cast<size_t>(file.dim1_lens);
                            size_t send_nelems = static_cast<size_t>(target_need) *
                                                 static_cast<size_t>(file.dim1_lens);

                            // 使用 world_rank 作为 tag，确保唯一性
                            int unique_tag = t_file_idx * file.world_procs_num + wr;

                            // if (target != rank_local)
                            // {
                            MPI_Request req;
                            MPI_Isend(data_block + send_elem_start,
                                      static_cast<int>(send_nelems),
                                      MPI_FLOAT,
                                      target,
                                      unique_tag, // tags
                                      color_comm, &req);
                            send_reqs.push_back(req);
                            // }
                        }
                    }
                }
            }

            // cout << "[INFO] Process " << world_rank << " posted all Isend requests for variable '" << varname << "'\n";
            // 等待所有通信完成
            if (!send_reqs.empty())
                MPI_Waitall(static_cast<int>(send_reqs.size()), send_reqs.data(), MPI_STATUSES_IGNORE);
            if (!recv_reqs.empty())
                MPI_Waitall(static_cast<int>(recv_reqs.size()), recv_reqs.data(), MPI_STATUSES_IGNORE);

            // cout << "[INFO] Process " << world_rank << " completed all communication for variable '" << varname << "'\n";
            // 排序并写入 NetCDF

            // if (rank_local == 0)
            // {
            //     // cout << "[INFO] Process " << world_rank << " starting to write NetCDF for variable '" << varname << "'\n";
            //     BinFileData &first_file = bin_files_data[0];
            //     for (int i = 0; i < first_file.mycount; ++i)
            //     {
            //         for(int j = 0; j < first_file.dim1_lens; ++j)
            //         {
            //             cout << recv_buffers[0][(size_t)i * (size_t)first_file.dim1_lens + (size_t)j] << " ";
            //         }
            //         cout << "\n";
            //     }
            // }

            int ncid, ierr;
            string s_outdir = outdir ? outdir : "";
            // string nc_file = s_outdir + "_" + varname + ".nc";
            string nc_file = "./"  + varname + ".dir/" +
                             s_outdir + "_" + varname + ".nc";

            ierr = ncmpi_open(color_comm, nc_file.c_str(), NC_WRITE, MPI_INFO_NULL, &ncid);
            if (ierr != NC_NOERR)
            {
                printf("nc_open_par error: %s\n", ncmpi_strerror(ierr));
            }
            int nc_varid;
            ierr = ncmpi_inq_varid(ncid, varname.c_str(), &nc_varid);

            for (int file_idx = 0; file_idx < FILE_NUM; ++file_idx)
            {
                if (file_idx < bin_files_data.size())
                {
                    BinFileData &file = bin_files_data[file_idx];
                    vector<float> put_buffer((size_t)file.mycount * (size_t)file.dim1_lens, 0.0f);
                    // reorder
                    for (int i = 0; i < file.mycount; ++i)
                    {
                        int lidx = file.idx_recv[i] - 1;
                        if (lidx < 0 || lidx >= file.mycount)
                        {
                            fprintf(stderr, "[ERROR] idx_recv[%d] = %d out of range [0, %d)\n", i, lidx, file.mycount);
                        }
                        else
                        {
                            const float *src = recv_buffers[file_idx].data() + (size_t)i * (size_t)file.dim1_lens;
                            float *dest = put_buffer.data() + (size_t)lidx * (size_t)file.dim1_lens;
                            memcpy(dest, src, (size_t)file.dim1_lens * sizeof(float));
                        }
                    }
                    MPI_Offset start[2] = {static_cast<MPI_Offset>(file.mystart - 1), 0};
                    MPI_Offset count[2] = {static_cast<MPI_Offset>(file.mycount), static_cast<MPI_Offset>(file.dim1_lens)};
                    ierr = ncmpi_put_vara_float_all(ncid, nc_varid, start, count, put_buffer.data());
                    if (ierr != NC_NOERR)
                    {
                        printf("PnetCDF error: %s\n", ncmpi_strerror(ierr));
                    }
                }
                else
                {
                    // 没有数据也要参与 collective 写入
                    MPI_Offset start[2] = {0, 0};
                    MPI_Offset count[2] = {0, 0};
                    ierr = ncmpi_put_vara_float_all(ncid, nc_varid, start, count, nullptr);
                    if (ierr != NC_NOERR)
                    {
                        printf("PnetCDF error (empty): %s\n", ncmpi_strerror(ierr));
                    }
                }
            }
            // 清理
            send_reqs.clear();
            recv_reqs.clear();
            for (auto &file : bin_files_data)
                for (auto ptr : file.data_ptrs)
                    delete[] ptr;
            bin_files_data.clear();
            ncmpi_close(ncid);
            MPI_Barrier(color_comm);
            cout << "[INFO] Process " << world_rank << " finished writing NetCDF for variable '" << varname << "'\n";
        }
    }
    MPI_Comm_free(&color_comm);
}

// 主函数
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int world_size = -1, world_rank = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (argc >= 2)
    {
        FILE_NUM = stoi(argv[1]);
        // cout << "[INFO] Using FILE_NUM = " << FILE_NUM << endl;
        if (FILE_NUM <= 0)
        {
            FILE_NUM = 16;
            cout << "[WARN] Invalid FILE_NUM specified. Using default value: " << FILE_NUM << endl;
        }
    }
    else
    {
        // 当未指定 FILE_NUM 时，保持默认值或后续根据文件数动态调整
        if (world_rank == 0)
        {
            cout << "[INFO] FILE_NUM not specified. Will be determined automatically." << endl;
        }
    }

    // printf("MPI initialized.\n");
    string info_path = "scalar_1d_list_info.txt";
    ifstream fin(info_path);
    if (!fin)
    {
        cerr << "Cannot open " << info_path << endl;
        // MPI_Finalize();
        return 1;
    }

    string line;
    while (getline(fin, line))
    {
        // 跳过空行
        if (line.empty())
            continue;
        // 检查是否是新的一组
        if (line.substr(0, 7) != "ncfile=")
            continue;

        string ncfile_str = line.substr(7); // 跳过'ncfile='

        // 提取时间步字段
        string timestep_dir;
        smatch m;
        regex re(R"((\d{4}-\d{2}-\d{2}-\d+))");
        if (regex_search(ncfile_str, m, re))
        {
            timestep_dir = m[1];
        }
        else
        {
            timestep_dir = "unknown_time";
        }
        // printf("[INFO] Processing ncfile: %s (time step: %s)\n", ncfile_str.c_str(), timestep_dir.c_str());

        // 读取变量个数
        if (!getline(fin, line))
            break;
        int nvars = stoi(line.substr(6)); // 跳过'nvars='

        // 跳过'varnames:'行
        if (!getline(fin, line))
            break;

        // 读取变量名
        vector<string> varnames;
        for (int i = 0; i < nvars; ++i)
        {
            if (!getline(fin, line))
                break;
            varnames.push_back(line);
        }

        // 构造C字符串
        vector<char> varnames_c;
        for (const auto &varname : varnames)
        {
            varnames_c.insert(varnames_c.end(), varname.begin(), varname.end());
            varnames_c.push_back('\0');
        }
        const char *varnames_c_ptr = varnames_c.data();

        if (world_rank == 0)
        {
            // 打印结果
            cout << "ncfile: " << ncfile_str << endl;
            cout << "nvars: " << nvars << endl;
            cout << "varnames:" << endl;
            for (const auto &v : varnames)
            {
                cout << "  " << v << endl;
            }
            cout << varnames_c_ptr << endl;

            // 检查varnames_c_ptr内容
            cout << "varnames_c_ptr逐个变量名内容:" << endl;
            const char *p = varnames_c_ptr;
            for (int i = 0; i < nvars; ++i)
            {
                cout << "  [" << i << "]: " << p << endl;
                p += strlen(p) + 1;
            }
        }

        // 取第一个变量名
        int bin_count = 0;
        if (!varnames.empty())
        {
            string first_var = varnames[0];
            string s_outdir = ncfile_str;
            for (int i = 0;; ++i)
            {
                // string bin_path = s_outdir + "_" + first_var + "_grp" + to_string(i) + ".bin";
                string bin_path = "./"  + first_var + ".dir/" +
                                  s_outdir + "_" + first_var + "_grp" +
                                  to_string(i) + ".bin";

                if (ifstream(bin_path).good())
                {
                    ++bin_count;
                }
                else
                {
                    break; // 连续编号，遇到不存在就停止
                }
            }
            // cout << "Bin file count for first variable: " << bin_count << endl;
        }

        // 如果未从命令行获取 FILE_NUM，则动态计算
        if (argc < 2 && bin_count > 0)
        {
            FILE_NUM = (bin_count + world_size - 1) / world_size;
            if (FILE_NUM == 0) FILE_NUM = 1; // 确保至少为1
            if (world_rank == 0)
            {
                cout << "[INFO] Auto-calculated FILE_NUM = " << FILE_NUM 
                     << " (bin_count=" << bin_count << ", world_size=" << world_size << ")" << endl;
            }
        }

        // 创建一个只包含前 bin_count 个进程的新 communicator
        MPI_Comm sub_comm = MPI_COMM_NULL;

        // 计算工作进程数
        int work_size = (bin_count + FILE_NUM - 1) / FILE_NUM;
        if (work_size > world_size)
        {
            if (world_rank == 0)
            {
                fprintf(stderr, "[ERROR] 需要进程数 %d，但实际只有 %d 个进程，无法处理所有 bin 文件！\n", work_size, world_size);
            }
            fin.close();
            MPI_Finalize();
            return 1;
        }

        int color = (world_rank < work_size) ? 0 : MPI_UNDEFINED;
        MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &sub_comm);
        if (color == 0)
        {
            printf("[INFO] Process %d in sub_comm handling ncfile: %s\n", world_rank, ncfile_str.c_str());
            // 只在新 comm 内的进程参与转换
            bin_to_nc_group_sp(timestep_dir, ncfile_str.c_str(), MPI_Comm_c2f(sub_comm), 0, nvars, varnames_c_ptr);
            MPI_Comm_free(&sub_comm);
        }
    }
    fin.close();
    MPI_Finalize();
    return 0;
}
