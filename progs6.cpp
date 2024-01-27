#include "mpi.h"
#include <iostream>
#include <chrono>

using namespace std;

int main ( int argc, char *argv[] ) {
    long N=10;
    if (argc > 1) {
        N = atoi(argv[1]);
        argc-=1;
    }

    int ProcNum, ProcRank, RecvRank;

    MPI_Status Status;
    MPI_Init ( &argc, &argv );

    MPI_Comm_size ( MPI_COMM_WORLD, &ProcNum);      //определяем количетсво процессов
    MPI_Comm_rank ( MPI_COMM_WORLD, &ProcRank);     // определяем текущий процесс

    if (ProcNum==1) {cout <<"ProcNum == 1, return" <<endl; return 0;}

    long k = ProcNum - 1; // поток 0 не участвует
    long n = (N + k - 1) / k;

    if (ProcRank == 0) {                            // 0 процесс (рассылка и принятие)
        std::srand(std::time(nullptr));

        N = n * k;
        u_int8_t byte_array[N];
        u_int8_t counter =0;
        for (long i=0; i<N; ++i){
            byte_array[i] = std::rand()%256;
        }

        std::chrono::time_point<std::chrono::system_clock> m_StartTime = std::chrono::system_clock::now();

        for (int i=0; i<10; ++i) {
            MPI_Sendrecv_replace(&byte_array[0], n, MPI_BYTE, 14, 0, 14, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        }
        
        std::chrono::time_point<std::chrono::system_clock> m_EndTime = std::chrono::system_clock::now();

        cout << N << "," << k << "," << std::chrono::duration_cast<std::chrono::microseconds>(m_EndTime - m_StartTime).count()/10 << ";";
    } else {                                        // 1 - ProcNum процессы (вычисления)
        if (ProcRank == 14)
        for (int i=0; i<10; ++i){
            u_int8_t my_enother_u_int8_t_array[n];
            //MPI_Sendrecv_replace() 
            MPI_Sendrecv_replace(&my_enother_u_int8_t_array[0], n, MPI_BYTE, 0, 0, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        }
    }

    MPI_Finalize();
    return 0;
}

