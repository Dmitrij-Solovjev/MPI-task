#include "mpi.h"
#include <iostream>
#include <chrono>

using namespace std;

const long N = 10; // 00000000;

u_int8_t byte_array[N];


u_int8_t search_min(long n, u_int8_t *my_enother_byte_array[]){
    u_int8_t min = 255;
    for (long i=0; i<n; ++i)
        if ((*my_enother_byte_array)[i] < min)
            min = (*my_enother_byte_array)[i];
    
    return min;
}

int main ( int argc, char *argv[] ) {
    int ProcNum, ProcRank, RecvRank;
    MPI_Status Status;
    MPI_Init ( &argc, &argv );
    MPI_Comm_size ( MPI_COMM_WORLD, &ProcNum);      //определяем количетсво процессов
    MPI_Comm_rank ( MPI_COMM_WORLD, &ProcRank);     // определяем текущий процесс

    if (ProcNum==1) {cout <<"ProcNum == 1, return" <<endl; return 0;}
    long k = ProcNum - 1;

    long n = (N + k - 1) / k;
    u_int8_t *my_enother_u_int8_t_array[n];
    if (ProcRank == 0) {                            // 0 процесс
        // разослать нужные участки массива.
        
        std::chrono::time_point<std::chrono::system_clock> m_StartTime = std::chrono::system_clock::now(), m_EndTime;

        for (int i=1; i<ProcNum; ++i){
            MPI_Send(&byte_array[i*n], n, MPI_BYTE, i, 0, MPI_COMM_WORLD);  
        }
    
        u_int8_t min = 255, check_min;
        for (int i=1; i<ProcNum; ++i){
            MPI_Recv(&check_min, 1, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            if (check_min < min) min=check_min;
        }

    } else {
        MPI_Recv(my_enother_u_int8_t_array[0], n, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
        u_int8_t min = search_min(n, &my_enother_u_int8_t_array[0]);
        MPI_Send(&min, 1, MPI_BYTE, 0, 0, MPI_COMM_WORLD);  
    }

    cout << ProcNum << " " << ProcRank <<  endl;
    MPI_Finalize();
    return 0;
}