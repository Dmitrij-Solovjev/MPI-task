#include "mpi.h"
#include <iostream>
#include <chrono>
#include <unistd.h>

using namespace std;

int main ( int argc, char *argv[] ) {
    long N=10;
    double sleep_time = 5;
    if (argc > 1){
        N = atoi(argv[1]);
        sleep_time = atof(argv[2]);
        argc-=2;
//        cout << N << " " << sleep_time << endl;
    }
    sleep_time/=100.0;

    int ProcNum, ProcRank, RecvRank;

    MPI_Status Status;
    MPI_Init ( &argc, &argv );

    MPI_Comm_size ( MPI_COMM_WORLD, &ProcNum);      //определяем количетсво процессов
    MPI_Comm_rank ( MPI_COMM_WORLD, &ProcRank);     // определяем текущий процесс

    if (ProcNum==1) {cout <<"ProcNum == 1, return" <<endl; return 0;}

    if (ProcRank == 0) {                            // 0 процесс (рассылка и принятие)
        u_int8_t byte_array[N];

        std::chrono::time_point<std::chrono::system_clock> m_StartTime = std::chrono::system_clock::now();

        // разослать нужные участки массива.
        for (int i=1; i<ProcNum; ++i) {
            MPI_Send(&byte_array[0], N, MPI_BYTE, i, 0, MPI_COMM_WORLD);
        }

        u_int8_t min = 255, check_min;

        for (int i=1; i<ProcNum; ++i) {
            MPI_Recv(&byte_array[0], N, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        }

        std::chrono::time_point<std::chrono::system_clock> m_EndTime = std::chrono::system_clock::now();

        cout << N << "," << ProcNum-1 << "," << sleep_time << "," << std::chrono::duration_cast<std::chrono::microseconds>(m_EndTime - m_StartTime).count() << ";";
    
    } else {                                        // 1 - ProcNum процессы (вычисления)
        u_int8_t my_enother_u_int8_t_array[N];
        MPI_Recv(&my_enother_u_int8_t_array[0], N, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        
        sleep(sleep_time/(ProcNum-1));

        MPI_Send(&my_enother_u_int8_t_array[0], N, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

