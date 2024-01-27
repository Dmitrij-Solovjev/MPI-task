#include "mpi.h"
#include <iostream>
#include <chrono>
#include <iterator>
#include <algorithm>

using namespace std;

int ProcNum, ProcRank, RecvRank;
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
    if (ProcRank == 0) {

        int counter = 0;
        for (int i=0; i<ProcNum and counter<N; ++i) {
            cout << "send to:" << i << endl;
            int target_num =  (N - i + ProcNum-1) / ProcNum;
            
            MPI_Send(&A[0][0], N*N, MPI_BYTE, (i%ProcNum)+1, 0, MPI_COMM_WORLD);
            
            MPI_Send(&B[counter][0], N*target_num, MPI_BYTE, i%ProcNum+1, 0, MPI_COMM_WORLD);
            counter+=target_num;
        }

        int target_num_max =  (N + ProcNum-1) / ProcNum;
        uint8_t buff[target_num_max][N];
        for (int i=0; i<ProcNum; ++i) {
            MPI_Recv(&buff[0][0], N*target_num_max, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            
            int target_num =  (N - (Status.MPI_SOURCE-1) + ProcNum-1) / ProcNum;
            
            int position = N/ProcNum*(Status.MPI_SOURCE-1);  // целые перед нами
            cout << position << endl;
            if (N%ProcNum >= Status.MPI_SOURCE) {            // мы внутри остатка
                position+=(Status.MPI_SOURCE-1);
            } else {
                position+= N%ProcNum;               // мы правее остатка
            }
            cout << "For " << Status.MPI_SOURCE << " pos: " <<  position << endl;

            //copy(&buff[0][0], &buff[target_num][N-1], &C[position][0]);
            for (int j=0; j<target_num; ++j) {
                for (int k=0; k<N; ++k) {
                    C[position+j][k] = buff[j][k];
                }
            }
        }
    }
/*********************************************************************************************************************/
/*********************************************************************************************************************/
    else {
        
        // определение сколько будет выслано строк
        int target_num =  (N - (ProcRank-1) + ProcNum-1) / ProcNum;
        int target_num_max =  (N + ProcNum-1) / ProcNum;

        // Получение А[][]
        if (target_num == 0) return; 
        MPI_Recv(&A[0][0], N*N, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        
        // Получение B[]
        MPI_Recv(&B[0][0], N*target_num, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        cout << "Принял B" << endl;
        
        // Вычисление C[]
        for (int i=0; i<target_num; ++i) {
            vectormult(N, A, B[i], C[i]);
        }

        // Отправка С[]
        MPI_Send(&C[0][0], N*target_num_max, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
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
    //    cout << "lenght_of_matrix:" << N <<endl;
    }

    MPI_Init ( &argc, &argv );
    MPI_Comm_size ( MPI_COMM_WORLD, &ProcNum);      //определяем количетсво процессов
    MPI_Comm_rank ( MPI_COMM_WORLD, &ProcRank);     // определяем текущий процесс

    ProcNum--;                                      // для расчетов используется на 1 меньше

    if (ProcNum> N){
        ProcNum = N;
    }
    
    uint8_t **A, **B, **C;
    A = alloc_2d_int(N,N);
    B = alloc_2d_int(N,N);
    C = alloc_2d_int(N,N);
    
    if (ProcRank == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = rand() % 10;
            }
        }
    
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                B[i][j] = rand() % 10;
            }
        }
    }
    std::chrono::time_point<std::chrono::system_clock> m_StartTime = std::chrono::system_clock::now();

    mult_lin(N, A, B, C);

    if (ProcRank == 0) {
        cout << "final cout C" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << (int)C[i][j] << " ";
        }
        cout << endl;
    }
    cout << "FINISH" << endl;
    }

    if (ProcRank == 0) {
        std::chrono::time_point<std::chrono::system_clock> m_EndTime = std::chrono::system_clock::now();
        cout << N << "," << ProcNum << "," << std::chrono::duration_cast<std::chrono::microseconds>(m_EndTime - m_StartTime).count()/10 << ";";
    }

    MPI_Finalize();
    return 0;
}