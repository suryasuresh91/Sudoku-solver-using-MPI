#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>
#include <time.h>
#include "list.h"

#define UNASSIGNED 0
#define UNCHANGEABLE -1

#define POS 0
#define VAL 1

#define TAG_HYP     1
#define TAG_EXIT    2
#define TAG_ASK_JOB 3

#define ROW(i) i/m_size
#define COL(i) i%m_size
#define BOX(row, col) r_size*(row/r_size)+col/r_size

#define BLOCK_LOW(rank, p, n) ((rank)*(n)/(p))
#define BLOCK_HIGH(rank, p, n) (BLOCK_LOW(rank+1,p,n)-1)

void init_masks(int* sudoku, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask);
void update_masks(int num, int row, int col, uint64_t *rows_mask, uint64_t *cols_mask, uint64_t *boxes_mask);
void rm_num_masks(int num, int row, int col, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask);
int solving_sudoku(int* sudoku, int* cp_sudoku, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, List* work, int last_pos);
int is_safe_num( uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, int row, int col, int num);
int exists_in( int index, uint64_t* mask, int num);
void delete_from(int* sudoku, int *cp_sudoku, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, int cell);

void send_ring(void *msg, int tag, int dest);
int* read_matrix(char *argv[]);
void print_sudoku(int *sudoku);
int new_mask( int size);
int solve(int *sudoku);
Item invalid_hyp(void);

int r_size, m_size, v_size, rank, p;

int main(int argc, char *argv[]){
      clock_t begin = clock();

    int* sudoku, result, total;

    if(argc == 2){
        sudoku = read_matrix(argv);

        MPI_Init (&argc, &argv);
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size (MPI_COMM_WORLD, &p);

        result = solve(sudoku);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&result, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if(!total && !rank)
            printf("No solution\n");
        else if(total && result){
            printf("\n Rank = %d \n", rank);
            printf("\n     SOLUTION: \n\n");
            print_sudoku(sudoku);
        }

        fflush(stdout);
        MPI_Finalize();

    }
    else
        printf("invalid input arguments.\n");

    free(sudoku);

        clock_t end = clock();
        double execution_time = (double)(end - begin)/CLOCKS_PER_SEC;
        printf("\n ****Rank = %d --- Execution time : %f microseconds\n", rank, execution_time);

    return 0;
}

int solve(int* sudoku){
    int i, last_pos, flag_start = 0, solved = 0;
    Item hyp;

    uint64_t *rows_mask = (uint64_t*) malloc(m_size * sizeof(uint64_t));
    uint64_t *cols_mask = (uint64_t*) malloc(m_size * sizeof(uint64_t));
    uint64_t *boxes_mask = (uint64_t*) malloc(m_size * sizeof(uint64_t));
    int *cp_sudoku = (int*) malloc(v_size * sizeof(int));
    //A work list contains the pairs (cell_id, number) from which serial DFS search must still be performed
    //(in other words it contains the root values of the unexplored parts of the search tree)
    List *work = init_list();

    for(i = 0; i < v_size; i++) {
        if(sudoku[i])
            cp_sudoku[i] = UNCHANGEABLE;
        else{
            cp_sudoku[i] = UNASSIGNED;
            if(!flag_start){
                flag_start = 1;
                hyp.cell = i;
            }
            last_pos = i;
        }
    }

    init_masks(sudoku, rows_mask, cols_mask, boxes_mask);

    //calculate the low and high values for the first cell for each process
    //and insert it in the work list
    //Assign a different set of numbers to each process to find their possible locations i.e. different work lists.
    for(i = 1 + BLOCK_HIGH(rank, p, m_size); i >= 1 + BLOCK_LOW(rank, p, m_size); i--){
        hyp.num = i;
        insert_head(work, hyp);
    }

    // try to solve sudoku
    solved = solving_sudoku(sudoku, cp_sudoku, rows_mask, cols_mask, boxes_mask, work, last_pos);

    if(solved){
        //if the solution is found copy the solution to be retrieved
        for(i = 0; i < v_size; i++)
            if(cp_sudoku[i] != UNCHANGEABLE)
                sudoku[i] = cp_sudoku[i];
    }

    free(work);
    free(rows_mask);
    free(cols_mask);
    free(boxes_mask);
    free(cp_sudoku);

    return solved;
}

int solving_sudoku(int* sudoku, int* cp_sudoku, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, List* work, int last_pos){
    int i, cell, val, number_amount, flag = 0, no_sol_count, len;

    MPI_Request request;
    MPI_Status status;
    Item hyp, no_hyp = invalid_hyp();

    //a while loop to get work from other processes when the list gets empty
    while(1){

        //a while loop to get work from the work list
        while(work->head != NULL){

            //pop a probable number from the work list
            hyp = pop_head(work);
            len = work->len;
            int start_pos = hyp.cell;

            //if the number of the initial hypothesis is not valid skip this hypothesis
            //check if it is safe to add number in hyp.cell i.e hyp.num already exists in row/col/box
            if(!is_safe_num(rows_mask, cols_mask, boxes_mask, ROW(hyp.cell), COL(hyp.cell), hyp.num))
                continue;

            //a while loop to solve the sudoku
            while(1){

                //listen to incoming messages
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

                //if a message has been received
                if(flag && status.MPI_TAG != -1){
                    flag = 0;

                    //find the size of the message and allocate a buffer for the message
                    MPI_Get_count(&status, MPI_INT, &number_amount);
                    int* number_buf = (int*)malloc(number_amount * sizeof(int));

                    //read the message to the allocated buffer
                    MPI_Recv(number_buf, number_amount, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                    //if the message is an exit signal, forward the message in the ring and return
                    if(status.MPI_TAG == TAG_EXIT){  //TAG_EXIT = 2
                        send_ring(&rank, TAG_EXIT, -1);
                        return 0;
                    }
                                        //if the message is a job request
                    else if(status.MPI_TAG == TAG_ASK_JOB){ //TAG_ASK_JOB = 3

                        //and if there is work to give in the list
                        if(work->tail != NULL){

                            //remove an hypothesis from the list
                            Item hyp_send = pop_tail(work);

                            //concatenate the hypothesis with the sudoku and send it
                            int* send_msg = (int*)malloc((v_size+2)*sizeof(int));
                            memcpy(send_msg, &hyp_send, sizeof(Item));
                            memcpy((send_msg+2), cp_sudoku, v_size*sizeof(int));

                            //send the message
                            MPI_Send(send_msg, (v_size+2), MPI_INT, status.MPI_SOURCE, TAG_HYP, MPI_COMM_WORLD);
                            free(send_msg);
                        }
                        else //if there isn't work to do send an impossible hypothesis message
                            MPI_Send(&no_hyp, 2, MPI_INT, status.MPI_SOURCE, TAG_HYP, MPI_COMM_WORLD);
                    }
                    free(number_buf);
                }

                //update the masks and sudoku with the hypothesis removed from the list
                update_masks(hyp.num, ROW(hyp.cell), COL(hyp.cell), rows_mask, cols_mask, boxes_mask);
                cp_sudoku[hyp.cell] = hyp.num;

                //iterate cells of the sudoku
                for(cell = hyp.cell + 1; cell < v_size; cell++){

                    //find a subsequent cell which does not have a value yet
                    if(cp_sudoku[cell]) //if the cell has an unchangeable number skip the cell
                        continue;

                    //test numbers in the cell////find all possibles value which can go in that cell
                    for(val = m_size; val >= 1; val--){

                        //if the current number is not valid in this cell skip the number
                        if(is_safe_num(rows_mask, cols_mask, boxes_mask, ROW(cell), COL(cell), val)){

                            //if the cell is the last one and a valid number for it was found the sudoku has been solved
                            //send an exit signal message to the ring of communication and return
                            if(cell == last_pos){
                                cp_sudoku[cell] = val;
                                send_ring(&rank, TAG_EXIT, -1);
                                return 1;
                            }

                            //insert the safe number for the cell as an hypothesis in the work list
                            hyp.cell = cell;
                            hyp.num = val;
                            insert_head(work, hyp);
                        }

                    }

                    break;
                }

                if(work->len == len){
                    for(cell = v_size - 1; cell >= start_pos; cell--)
                        if(cp_sudoku[cell] > 0){
                            rm_num_masks(cp_sudoku[cell],  ROW(cell), COL(cell), rows_mask, cols_mask, boxes_mask);
                            cp_sudoku[cell] = UNASSIGNED;
                        }
                    break;
                }

                //take a new hypothesis from the work list
                hyp = pop_head(work);

                //clear the sudoku down to the point of that hypothesis
                for(cell--; cell >= hyp.cell; cell--){
                    if(cp_sudoku[cell] > 0) {
                        rm_num_masks(cp_sudoku[cell],  ROW(cell), COL(cell), rows_mask, cols_mask, boxes_mask);
                        cp_sudoku[cell] = UNASSIGNED;
                    }
                }
            }
        }

        no_sol_count = 0;

        if(p == 1)
            return 0;

        //cycle to ask other processes for work
        for(i = rank + 1;; i++){

            if(i == p) i = 0;
            if(i == rank) continue;

            //send a work request message to the ith process
            MPI_Send(&i, 1, MPI_INT, i, TAG_ASK_JOB, MPI_COMM_WORLD);

            //wait for an incoming message from any process
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            //find the size of the received message and allocate a buffer for the message
            MPI_Get_count(&status, MPI_INT, &number_amount);
            int* number_buf = (int*)malloc(number_amount * sizeof(int));

            //read the message to the allocated buffer
            MPI_Recv(number_buf, number_amount, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            //if the message is a new hypothesis
            if(status.MPI_TAG == TAG_HYP && number_amount != 2){

                //take the hypothesis and sudoku information from the message
                Item hyp_recv;
                memcpy(&hyp_recv, number_buf, sizeof(Item));
                memcpy(cp_sudoku, (number_buf+2), v_size*sizeof(int));

                //delete everything to the point of the received hypothesis
                delete_from(sudoku, cp_sudoku, rows_mask, cols_mask, boxes_mask, hyp_recv.cell);

                //insert it in the work list
                insert_head(work, hyp_recv);
                free(number_buf);

                break;

            //if the message is an invalid hypothesis increase the number of processos without work
            }else if(status.MPI_TAG == TAG_HYP && number_amount == 2){
                no_sol_count++;

            //if the message is an exit signal forward the signal to the ring and return
            }else if(status.MPI_TAG == TAG_EXIT){
                send_ring(&rank, TAG_EXIT, -1);
                return 0;

            //if the message is a request for work send a no work to give message
            }else if(status.MPI_TAG == TAG_ASK_JOB){
                MPI_Send(&no_hyp, 2, MPI_INT, status.MPI_SOURCE, TAG_HYP, MPI_COMM_WORLD);
            }

            //if all other processes had no work to give there is no solution to the sudoku
            if(no_sol_count == p-1 && rank == 0){
                send_ring(&rank, TAG_EXIT, -1);
                return 0;
            }

            free(number_buf);
        }
    }
}

//when a process receives a new hypothesis and the corresponding cp_sudoku, then we has to clear all the work
//that the sender process have done after the hypothesis's position
void delete_from(int* sudoku, int *cp_sudoku, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, int cell){
    int i;

    //initialize the masks from the sudoku read from the file
    init_masks(sudoku, rows_mask, cols_mask, boxes_mask);

    //clear all work done by the sender process
    for(i = v_size; i >= cell; i--)
        if(cp_sudoku[i] > 0)
            cp_sudoku[i] = UNASSIGNED;

    //initialize the masks from the cp_sudoku received
    for(i = 0; i < cell; i++)
        if(cp_sudoku[i] > 0)
            update_masks(cp_sudoku[i], ROW(i), COL(i), rows_mask, cols_mask, boxes_mask);
}

//create an invalid hypothesis
Item invalid_hyp(void){
    Item item;
    item.cell = -1;
    item.num = -1;
    return item;
}

void send_ring(void *msg, int tag, int dest){
    int msg_send[2];
    msg_send[0] =*((int*) msg);
    msg_send[1] = dest;

    if(rank == p-1)
        MPI_Send(msg_send, 2, MPI_INT, 0, tag, MPI_COMM_WORLD);
    else
        MPI_Send(msg_send, 2, MPI_INT, rank + 1, tag, MPI_COMM_WORLD);
}

//test if a number (num) exists in a row, column or box. To check if it exists in a row, the mask of that row is passed as input
int exists_in(int index, uint64_t* mask, int num) {
    int res, masked_num = 1 << (num-1);

    res = mask[index] | masked_num;  //example row_mask = 1010 ; number to check = 0100 (3) then row_masl | num = 1110
    if(res != mask[index])           // since the result is different from the initial mask, it is safe to put 3 in this row/col/box
        return 0;
    return 1;
}

//check if a given number is safe in a given cell (row,col) for that check safety in corresponding row, column and box
int is_safe_num(uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, int row, int col, int num) {
    return !exists_in(row, rows_mask, num) && !exists_in(col, cols_mask, num) && !exists_in(BOX(row, col), boxes_mask, num);
}

//initialize a mask
/*int new_mask(int size) {
    return (0 << (size-1));
}*/


//initialize the masks from the sudoku read from the file
void init_masks(int* sudoku, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask) {
    int i;

    for(i = 0; i < m_size; i++){
        rows_mask[i]  = UNASSIGNED;
        cols_mask[i]  = UNASSIGNED;
        boxes_mask[i] = UNASSIGNED;
    }

    for(i = 0; i < v_size; i++)
        if(sudoku[i])
            update_masks(sudoku[i], ROW(i), COL(i), rows_mask, cols_mask, boxes_mask);
}

//remove number from the masks
void rm_num_masks(int num, int row, int col, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask) {
    int num_mask = 1 << (num-1);
    rows_mask[row] ^= num_mask;
    cols_mask[col] ^= num_mask;
    boxes_mask[BOX(row, col)] ^= num_mask;
}

//add new number to the masks
void update_masks(int num, int row, int col, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask) {
    int new_mask = 1 << (num-1);  //convert number found to mask ex: if dim=4x4, 3 = 0010
    rows_mask[row] |= new_mask;      //to add the new number to the current row's mask use bitwise OR
    cols_mask[col] |= new_mask;      //ex row_mask = 0101 ; number to add: 0010 --> row_mask OR num = 0111
    boxes_mask[BOX(row, col)] |= new_mask;
}

//read the input file
int* read_matrix(char *argv[]) {
    FILE *fp;
    size_t characters, len = 1;
    char *line = NULL, aux[3];
    int i, j, k, l;

    if((fp = fopen(argv[1], "r+")) == NULL) {
        fprintf(stderr, "unable to open file %s\n", argv[1]);
        exit(1);
    }

    getline(&line, &len, fp);
    r_size = atoi(line);
    m_size = r_size *r_size;
    v_size = m_size * m_size;

    int* sudoku = (int*)malloc(v_size * sizeof(int));

    k = 0, l = 0;
    len = m_size * 2;
    for(i = 0; (characters = getline(&line, &len, fp)) != -1; i++){
        for (j = 0; j < characters; j++) {
            if(isdigit(line[j])){
                aux[l++] = line[j];
            }else if(l > 0){
                aux[l] = '\0';
                l = 0;
                sudoku[k++] = atoi(aux);
                memset(aux, 0, sizeof aux);
            }
        }
    }

    free(line);
    fclose(fp);

    return sudoku;
}

void print_sudoku(int *sudoku) {
    int i;

    for (i = 0; i < v_size; i++) {
        if(i%m_size != m_size - 1){
            printf("%2d ", sudoku[i]);
            if (i% r_size == r_size -1)
                printf("  |  ");
        }
        else{
            printf("%2d\n\n", sudoku[i]);
            //printf("\n");
            if (i%(m_size*r_size)==(m_size*r_size)-1){
              for (int j = 0; j<m_size; j++)
                  printf("----");
              printf("\n");
            }

        }
    }
}
