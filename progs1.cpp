#include "mpi.h"
#include <iostream>
#include <chrono>

using namespace std;
//long N = 1000000;

u_int8_t search_min(long n, u_int8_t my_enother_byte_array[]) {
    u_int8_t min = 255;
    for (long i=0; i<n; ++i)
        if (my_enother_byte_array[i] < min)
            min = my_enother_byte_array[i];
    return min;
}

int main ( int argc, char *argv[] ) {
    long N=10;
    if (argc > 1){
        N = atoi(argv[1]);
//      cout << "N: " << N <<endl;
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
//        cout << "First proc! All=" << ProcNum << endl;

        N = n * k;
        u_int8_t byte_array[N];
        u_int8_t counter =0;
        for (long i=0; i<N; ++i){
            byte_array[i]=255;
            if (i % n == n - 1)
                byte_array[i] = counter++;
//            cout << (int)byte_array[i] << " " << endl;
        }
//        cout << "First proc! All=" << ProcNum << endl;

        std::chrono::time_point<std::chrono::system_clock> m_StartTime = std::chrono::system_clock::now();

        // разослать нужные участки массива.
        for (int i=1; i<ProcNum; ++i) {
            MPI_Send(&byte_array[(i-1)*n], n, MPI_BYTE, i, 0, MPI_COMM_WORLD);
        }
//        cout << "First proc! All=" << ProcNum << endl;


        u_int8_t min = 255, check_min;
        for (int i=1; i<ProcNum; ++i) {
            MPI_Recv(&check_min, 1, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
//            cout << (int)check_min << endl;
            if (check_min < min) min=check_min;
        }

        std::chrono::time_point<std::chrono::system_clock> m_EndTime = std::chrono::system_clock::now();
//        cout<<"v:" << (int)min << " t:" << std::chrono::duration_cast<std::chrono::microseconds>(m_EndTime - m_StartTime).count() << endl;
        cout << N << " " << k << " " << std::chrono::duration_cast<std::chrono::microseconds>(m_EndTime - m_StartTime).count() << endl;
    } else {                                        // 1 - ProcNum процессы (вычисления)
        u_int8_t my_enother_u_int8_t_array[n];
        MPI_Recv(&my_enother_u_int8_t_array[0], n, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        u_int8_t min = search_min(n, &my_enother_u_int8_t_array[0]);
        MPI_Send(&min, 1, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }

//    cout << ProcNum << " " << ProcRank <<  endl;
    MPI_Finalize();
    return 0;
}

