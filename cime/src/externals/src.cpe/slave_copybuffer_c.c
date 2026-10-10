#include <crts.h>

typedef struct{
    double *a, *b;
    int length, length1;
} str_par;

__thread_local crts_rply_t get_reply, put_reply;

void slave_copybuffer_(str_par* p)
{
    /* get pointer from fortran */
    double *ptr_a = p->a;
    double *ptr_b = p->b;
    /* get threads id */
    int id = CRTS_rid*8 + CRTS_cid;
    /* ldm buffer */
    int N = p->length;
    int i, jump;
    jump = N;
    if(id == 63)
    {
        N= p->length + p->length1;
    }
    double a[N], b[N];
    get_reply = 0;
    CRTS_dma_iget(a, ptr_a + id*jump, sizeof(double)*N, &get_reply);
    CRTS_dma_wait_value(&get_reply, 1);
    for(i = 0; i < N; i++)
    {
        b[i] = a[i];
    }
    put_reply = 0;
    CRTS_dma_iput(ptr_b + id*jump, b, sizeof(double)*N, &put_reply);
    CRTS_dma_wait_value(&put_reply, 1);
}