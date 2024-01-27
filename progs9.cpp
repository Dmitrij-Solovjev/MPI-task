#include "mpi.h"
#include <iostream>
#include <chrono>
#include <iterator>
#include <algorithm>

using namespace std;

int ProcNum, ProcRank, RecvRank, GigaProcNum;
MPI_Status Status;

void vectormult(int N, uint8_t *A[], uint8_t *b, uint8_t *c){
    for (int i = 0; i < N; ++i) {
        c[i] = 0;
        for (int j = 0; j < N; ++j) {
            c[i] += A[j][i] * b[j];
        }
    }
}

void mult_lin(int N, uint8_t *A[], uint8_t *B[], uint8_t *C[]) {
    int scatter_counts[ProcNum+1];                // количество
    int scatter_displs[ProcNum+1];                // смещение
    int counter = 0;
    scatter_counts[0]=0;
    scatter_displs[0]=0;
    
    for (int i=0; i<ProcNum; ++i) {
        int target_num =  (N - i + ProcNum-1) / ProcNum;
        scatter_counts[i+1]=target_num*N;
        scatter_displs[i+1]=counter*N;
        counter+=target_num;
    }
    for (int i=ProcNum+1; i<GigaProcNum; ++i){      // обрабатываем закрытые потоки (строк меньше чем потоков)
        scatter_counts[i]=0;
        scatter_displs[i]=0;
    }

    if (ProcRank == 0) {
        MPI_Bcast(&A[0][0], N*N, MPI_BYTE, 0, MPI_COMM_WORLD );
        MPI_Scatterv(
                    &B[0][0],
                    &scatter_counts[0],
                    &scatter_displs[0],
                    MPI_UINT8_T,
                    &B[0][0],
                    0,
                    MPI_UINT8_T,
                    0,
                    MPI_COMM_WORLD);

        MPI_Gatherv(
                    &C[0][0],
                    0,
                    MPI_UINT8_T,
                    &C[0][0],
                    &scatter_counts[0],
                    &scatter_displs[0],
                    MPI_UINT8_T,
                    0,
                    MPI_COMM_WORLD
        );
    }
/*********************************************************************************************************************/
/*********************************************************************************************************************/
    else {
        // определение сколько будет выслано строк
        int target_num =  (N - (ProcRank-1) + ProcNum-1) / ProcNum;

        
        if (target_num == 0 or ProcRank > ProcNum) return; 
        // Получение А[][]

        MPI_Bcast(&A[0][0], N*N, MPI_BYTE, 0, MPI_COMM_WORLD );

        MPI_Scatterv(
                    &B[0][0],
                    &scatter_counts[0],
                    &scatter_displs[0],
                    MPI_UINT8_T,
                    &B[0][0],
                    N*target_num,
                    MPI_UINT8_T,
                    0,
                    MPI_COMM_WORLD);
        
        // Вычисление C[]
        for (int i=0; i<target_num; ++i) {
            vectormult(N, A, B[i], C[i]);
        }

        MPI_Gatherv(
                    &C[0][0],
                    N*target_num,
                    MPI_UINT8_T,
                    &C[0][0],
                    &scatter_counts[0],
                    &scatter_displs[0],
                    MPI_UINT8_T,
                    0,
                    MPI_COMM_WORLD
        );
    }
}

uint8_t **alloc_2d_int(int rows, int cols) {
    uint8_t *data = (uint8_t *)malloc(rows*cols*sizeof(uint8_t));
    uint8_t **array= (uint8_t **)malloc(rows*sizeof(uint8_t*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

int main ( int argc, char *argv[] ) {
    long N=10;
    if (argc > 1) {
        N = atoi(argv[1]);
        argc-=1;
    }

    MPI_Init ( &argc, &argv );
    MPI_Comm_size ( MPI_COMM_WORLD, &GigaProcNum);      //определяем количетсво процессов
    MPI_Comm_rank ( MPI_COMM_WORLD, &ProcRank);     // определяем текущий процесс

    ProcNum = GigaProcNum - 1;                                      // для расчетов используется на 1 меньше

    if (ProcNum> N) {
        ProcNum = N+1;
    }
    
    uint8_t **A, **B, **C;
    A = alloc_2d_int(N,N);
    B = alloc_2d_int(N,N);
    C = alloc_2d_int(N,N);
    
    if (ProcRank == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = rand() % 10;
    //            cout << (int)A[i][j] << " ";
            }
    //        cout << endl; 
        }
    
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                B[i][j] = rand() % 10;
    //            cout << (int)B[i][j] << " ";
            }
    //        cout << endl;
        }
    }
    std::chrono::time_point<std::chrono::system_clock> m_StartTime = std::chrono::system_clock::now();

    mult_lin(N, A, B, C);

    if (ProcRank == 0) {
        std::chrono::time_point<std::chrono::system_clock> m_EndTime = std::chrono::system_clock::now();
        cout << N << "," << ProcNum << "," << std::chrono::duration_cast<std::chrono::microseconds>(m_EndTime - m_StartTime).count()/10 << ";";
    }

    MPI_Finalize();
    return 0;
}