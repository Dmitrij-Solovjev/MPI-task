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

uint64_t multiply(uint64_t n, u_int8_t my_enother_byte_array[]){
    uint64_t sum;
    for (uint64_t i=0; i<n; ++i)
        sum += my_enother_byte_array[i]* my_enother_byte_array[i+n];
    return sum;
}

int main ( int argc, char *argv[] ) {
    long N=10;
    if (argc > 1) {
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
        std::srand(std::time(nullptr));

        N = n * k * 2; // т.к. вектор
        u_int8_t byte_array[N];
        u_int8_t counter =0;
        for (long i=0; i<N; ++i){
            byte_array[i] = std::rand()%256;
        }
//        cout << "First proc! All=" << ProcNum << endl;

        std::chrono::time_point<std::chrono::system_clock> m_StartTime = std::chrono::system_clock::now();

        // разослать нужные участки массива.
        for (int i=1; i<ProcNum; ++i) {
            MPI_Send(&byte_array[(i-1)*2*n], 2*n, MPI_BYTE, i, 0, MPI_COMM_WORLD);
        }
//        cout << "First proc! All=" << ProcNum << endl;


        uint64_t sum = 0, gotten_value;
        for (int i=1; i<ProcNum; ++i) {
            MPI_Recv(&gotten_value, 1, MPI_UINT64_T, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            sum+=gotten_value;
            //if (check_min < min) min=check_min;
        }

        cout << sum << endl;
        std::chrono::time_point<std::chrono::system_clock> m_EndTime = std::chrono::system_clock::now();
//        cout<<"v:" << (int)min << " t:" << std::chrono::duration_cast<std::chrono::microseconds>(m_EndTime - m_StartTime).count() << endl;
        cout << N << " " << k << " " << std::chrono::duration_cast<std::chrono::microseconds>(m_EndTime - m_StartTime).count() << endl;
    } else {                                        // 1 - ProcNum процессы (вычисления)
        u_int8_t my_enother_u_int8_t_array[2*n];
        MPI_Recv(&my_enother_u_int8_t_array[0], 2*n, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        u_int8_t min = multiply(n, &my_enother_u_int8_t_array[0]);
        MPI_Send(&min, 1, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }

//    cout << ProcNum << " " << ProcRank <<  endl;
    MPI_Finalize();
    return 0;
}

