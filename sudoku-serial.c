#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "list.h"

#define UNASSIGNED 0
#define UNCHANGEABLE -1
#define ROW(i) i/m_size
#define COL(i) i%m_size
#define BOX(row, col) r_size*(row/r_size)+col/r_size

int solving_sudoku(int* sudoku, int* cp_sudoku, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, List* work, int last_pos);
void update_masks(int num, int row, int col, uint64_t *rows_mask, uint64_t *cols_mask, uint64_t *boxes_mask);
void rm_num_masks(int num, int row, int col, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask);
int is_safe_num( uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, int row, int col, int num);
int exists_in( int index, uint64_t* mask, int num);
int* read_matrix(char *argv[]);
void print_sudoku(int *sudoku);
int new_mask( int size);
int solve(int *sudoku);

int r_size, m_size, v_size;
int nr_iterations=0;

int main(int argc, char *argv[]){

    clock_t begin = clock();

    int* sudoku;

    if(argc == 2){
        sudoku = read_matrix(argv);
        printf("\n     PROBLEM : \n\n");
        print_sudoku(sudoku);

        if(solve(sudoku)){
              printf("\n     SOLUTION: \n\n");
              print_sudoku(sudoku);
        }else
            printf("No solution\n");
    }else
        printf("invalid input arguments.\n");

    free(sudoku);

    clock_t end = clock();
    double execution_time = (double)(end - begin)/CLOCKS_PER_SEC;
    printf("\n ****Execution time : %f microseconds\n\n", execution_time);

    return 0;
}

int solve(int* sudoku){
    int i, last_pos, flag_start = 0, solved = 0;
    Item hyp;

    uint64_t *rows_mask = (uint64_t*) malloc(m_size * sizeof(uint64_t));
    uint64_t *cols_mask = (uint64_t*) malloc(m_size * sizeof(uint64_t));
    uint64_t *boxes_mask = (uint64_t*) malloc(m_size * sizeof(uint64_t));
    int *cp_sudoku = (int*) malloc(v_size * sizeof(int));
    //A work list contains the pairs (cell_id, number) from wich serial DFS search is performed
    //(in other words it cointains the root values of the unexplored parts of the search tree)

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

    for(i = 0; i < m_size; i++){
        rows_mask[i]  = UNASSIGNED;
        cols_mask[i]  = UNASSIGNED;
        boxes_mask[i] = UNASSIGNED;
    }
//Initiate rows_mask, cols_mask, boxes_mask with sudoko values
    for(i = 0; i < v_size; i++)
        if(sudoku[i])
            update_masks(sudoku[i], ROW(i), COL(i), rows_mask, cols_mask, boxes_mask);
//insert all possible numbers into work list stack
    for(i = m_size; i >= 1; i--){
        hyp.num = i;
        insert_head(work, hyp);
    }

    solved = solving_sudoku(sudoku, cp_sudoku, rows_mask, cols_mask, boxes_mask, work, last_pos);
    if(solved) //not zero update sudoku with values from cp_sudoku
        for(i = 0; i < v_size; i++)
            if(cp_sudoku[i] != UNCHANGEABLE)
                sudoku[i] = cp_sudoku[i];

    free(work);
    free(rows_mask);
    free(cols_mask);
    free(boxes_mask);
    free(cp_sudoku);

    if(solved == 1)
        return 1;
    return 0;
}

int solving_sudoku(int* sudoku, int* cp_sudoku, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, List* work, int last_pos){
    int cell, val, len;
    Item hyp;

    //a while loop to get work from the work list
    while(work->head != NULL){
        hyp = pop_head(work); //pop a probable number from the work list
        len = work->len;
        int start_pos = hyp.cell;
        //check if it is safe to add number in hyp.cell i.e hyp.num already exists in row/col/box
        if(!is_safe_num(rows_mask, cols_mask, boxes_mask, ROW(hyp.cell), COL(hyp.cell), hyp.num))
            continue;
        //Once a safe number is found from work List, update masks& cp_sudoku
        //Else continue while. i.e., pop the next item from the stack
        while(1){
            nr_iterations++;
            //update the masks and sudoku with the hypothesis removed from the list
            update_masks(hyp.num, ROW(hyp.cell), COL(hyp.cell), rows_mask, cols_mask, boxes_mask);
            cp_sudoku[hyp.cell] = hyp.num;

            for(cell = hyp.cell + 1; cell < v_size; cell++){ //iterate cells of the sudoku
                //find a subsequent cell which does not have a value yet
                if(cp_sudoku[cell]) //if the cell has an unchangeable number skip the cell
                    continue;

                //find all possibles value which can go in that cell
                for(val = m_size; val >= 1; val--){
                    //if the current number is not valid in this cell skip the number
                    if(is_safe_num(rows_mask, cols_mask, boxes_mask, ROW(cell), COL(cell), val)){
                        //if the cell is the last one and a valid number for it was found the sudoku has been solved
                        if(cell == last_pos){
                            cp_sudoku[cell] = val;
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
                if(cp_sudoku[cell] > 0) { //remove interim values
                    rm_num_masks(cp_sudoku[cell],  ROW(cell), COL(cell), rows_mask, cols_mask, boxes_mask);
                    cp_sudoku[cell] = UNASSIGNED;
                }
            }
        }
    }
    return 0;
}

int exists_in(int index, uint64_t* mask, int num) {
    int res, masked_num = 1 << (num-1);  //int to mask of num

    res = mask[index] | masked_num;
    if(res != mask[index])
        return 0;
    return 1; //number already exists
}

int is_safe_num(uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask, int row, int col, int num) {
    return !exists_in(row, rows_mask, num) && !exists_in(col, cols_mask, num) && !exists_in(BOX(row, col), boxes_mask, num);
}

int new_mask(int size) {
    return (0 << (size-1));
}

void rm_num_masks(int num, int row, int col, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask) {
    int num_mask = 1 << (num-1);
    rows_mask[row] ^= num_mask;
    cols_mask[col] ^= num_mask;
    boxes_mask[BOX(row, col)] ^= num_mask;
}

void update_masks(int num, int row, int col, uint64_t* rows_mask, uint64_t* cols_mask, uint64_t* boxes_mask) {
    int new_mask = 1 << (num-1);
    rows_mask[row] |= new_mask;
    cols_mask[col] |= new_mask;
    boxes_mask[BOX(row, col)] |= new_mask;
}

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
