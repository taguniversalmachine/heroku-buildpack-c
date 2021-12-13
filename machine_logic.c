/*                                 
 * Machine Logic
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2013, TAG Universal Machine.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

 /* Random means YOU do the choosing.
 *  -- Ross Ashby
 *
 *
 *                                        
 *                                                        Cellular Automata Rule 30                                    Universal Constant
 *                                                                 Field of Action                                                           R30
 */


#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include "command.h"
#include "fmt_vec.h"
#include "cell_logic.h"
#include "machine_logic.h"
// Check the nth bit from the right end
#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
// CHECK_BIT(temp, n - 1)

/* Time Seed
 * A catalyst. 96 bits of environmental entropy in the form of
 * internal (Relative, Cyclic) and external (Absolute, Passing) time. 
 * Entropy is amplified to 128 bits to create the PSI fingerprint.
 * ts_cycle is the internal cycle counter to guarantee differentiation from call to call,
 * even if RTC clock does not increment.
 */

static uint32_t ts_cycle = 0;

TimeSeed time_seed()
{
  struct timeval tv; 
  ts_cycle += 1;
  if (gettimeofday(&tv, NULL) != 0) { halt("TimeSeed: gettimeofday failed");}               

  TimeSeed r = {.long_count = tv.tv_sec,    // RTC seconds
                .short_count = tv.tv_usec,  // RTC microseconds
                .cyclic = ts_cycle};        // Internal cyclic state

  return (r);
}

char* timeseed_tostr(TimeSeed* ts) {
   size_t  strsize = 100;
   char* r = malloc(strsize);
   memset(r, '\0', strsize);
   sprintf(r,"LC:%08X SC:%08X CY:%08X",ts->long_count, ts->short_count, ts->cyclic);
   return (r);
}

vec128bec_t* time_seed_to_vec(TimeSeed seed) {
   vec128bec_t* ts = vec_alloc();
   /* Convert 
      long count : long int (32 bits)
      short count: long int (32 bits)
      cyclic : long int     (32 bits)
      96 bits on the vector
   */

   int index = 1;
   // Long count
   for(int i=0; i<32; i++) {
     set_cell(ts, index + i, ((seed.long_count >> i) & 0x01) ? CELL_TRUE : CELL_FALSE);
   }

   // short count
   for(int i=0; i<32; i++){
    set_cell(ts, index + i + 32, ((seed.short_count >>i) & 0x01) ? CELL_TRUE : CELL_FALSE);
   }

   // cyclic
   for(int i=0; i<32; i++){
    set_cell(ts, index + i + 64, ((seed.cyclic >> 1) & 0x01) ? CELL_TRUE : CELL_FALSE);
   }

   // set remaining bits to CELL_FALSE
   for(int rest=96; rest<=128; rest++){
    set_cell(ts, rest, CELL_FALSE);
   }

   // every bit should be either CELL_TRUE or CELL_FALSE

   return (ts);
}

/* 
   Convert a string in PSI format [<:838087396B4405BCF017731EF1F99653:>] to vector
*/
int seed_str_to_vec(char* seed_str, vec128bec_t* out_vec) {
  int r = 0;
  char byte_string[3];
 /* printf("SEED_STR_TO_VEC\n"); */
  char* register_str;
  register_str = fmt_vecbe(out_vec, FMT_VEC_PSI);
 /* printf("%s\n", register_str); */
  if (seed_str == NULL) {
    vsetN(out_vec);
    return(-1);
  }

  if ((strlen(seed_str) == 38) && (seed_str[0] == '[') && (seed_str[1] == '<')  && (seed_str[2] == ':') &&
     (seed_str[35] == ':') && (seed_str[36] == '>')  && (seed_str[37] == ']'))
  {
    /* printf("SEED_STR_TO_VEC Looks Ok, converting bytes\n"); */
    for (int byte_num=0; byte_num < 16; byte_num++) {
      sprintf(byte_string, "%c%c",seed_str[(byte_num*2)+3], seed_str[((byte_num*2)+1)+3]);
   /*   printf("Byte:(%d):%s\n", byte_num, byte_string);  */
      int byte_int = (int)strtol(byte_string, NULL, 16);
  /*   printf("Converted to %d\n", byte_int);  */
      for (int bit=7; bit>=0; bit--){
        cell cell_val = CELL_FALSE;
        
        if (CHECK_BIT(byte_int, bit)) {
       /*   printf("CHECK_BIT %d,%d TRUE\n",byte_int, bit); */
          cell_val = CELL_TRUE;
        } else {
       /*   printf("CHECK_BIT %d,%d FALSE\n",byte_int, bit); */
          cell_val = CELL_FALSE;
        }
        int cell_num = 129 - ((byte_num * 8) + (7-bit) + 1);
     /*   printf("Setting cell %d to %d\n", cell_num, cell_val);  */
        set_cell(out_vec, cell_num, cell_val);
      }
    }
    return(0);
  } else {
    printf("SEED_STR_TO_VEC FORMAT CHECK FAILED, seed length was %ld\n", strlen(seed_str));
    vsetN(out_vec);
    return(-1);
  }
  return(r);
}

/* Compare TimeSeeds, return CELL_TRUE if equal, CELL_FALSE otherwise */
int ts_cmpE(TimeSeed* s0,  TimeSeed* s1){
  return ((s0->long_count  == s1->long_count)  &&
          (s0->short_count == s1->short_count) &&
          (s0->cyclic      == s1->cyclic)) ? CELL_TRUE : CELL_FALSE; 
}

void cp_setstate(struct cell_proc_t *restrict cp, cp_state_t st){
  cp->s = st;
}

cp_state_t cp_getstate(struct cell_proc_t *restrict cp) {
  return (cp->s);
}

void cp_halt(struct cell_proc_t *restrict cp, char* msg){
  cp_setstate(cp, CP_HALT);
  printf("Program Halt - %s.\n", msg);
  exit(EXIT_FAILURE);
}

/* Nullify memory slots and reset state
   C - does not deallocate any memory
 */
void cp_reset(struct cell_proc_t *restrict cp){
  cp->index = 0;
  cp->counter_mode = 0;
  cp->NR = cp->B;
  cp->CR = cp->A;

  cp_setstate(cp, CP_IDLE);
  frame_clear(cp->Stack);
  vsetN(cp->A);
  vsetN(cp->B);
  vsetN(cp->C);
  vsetN(cp->D);
  vsetN(cp->PSI);
  vsetN(cp->R30);
  vsetN(cp->R);
  vsetZ(cp->SDR30);
  vsetN(cp->SDTIME);
}

/* Allocate memory slots for registers and stack */
uint16_t cp_init(struct cell_proc_t *restrict cp){
  cp->A      = vec_alloc();
  cp->B      = vec_alloc();
  cp->C      = vec_alloc();
  cp->D      = vec_alloc();
  cp->PSI    = vec_alloc();
  cp->R30    = vec_alloc();
  cp->SDR30  = vec_alloc();
  cp->SDTIME = vec_alloc();
  cp->R      = vec_alloc();
  cp->Stack  = frame_alloc();
  cp_setstate(cp, CP_NULL);

  if (!((cp->A) && (cp->B) && (cp->C) && (cp->PSI) && (cp->R30) && (cp->R) )) {
     cp_halt(cp, "cp_init: Out of Memory");
  }

  cp_reset(cp);

  return(0);
}


/* Deallocate memory slots */
void cp_free(struct cell_proc_t *restrict cp){
  free(cp->A);
  free(cp->B);  
  free(cp->C);   
  free(cp->D);  
//  free(cp->X);
  free(cp->PSI); 
  free(cp->R30); 
  free(cp->SDR30); 
  free(cp->SDTIME); 
  free(cp->R);   
  frame_free(cp->Stack);
}

/* Machine Instructions */

/* CMPZ - Return CELL_TRUE if all cells are False, CELL_FALSE otherwise */
cell_state_t mi0_cmpZ (vec128bec_t* v) {
  uint8_t i;
  cell c;
  cell_state_t r = CELL_NULL;

  for(i = 1; i <= 128; i++) {
     c = get_cell(v, i);
     if ((c == CELL_TRUE) || (c == CELL_NULL)) {
        return (CELL_FALSE);
     }
  }

  r = CELL_TRUE;

  return r;
}

/* CMPN - Return CELL_TRUE if ANY cell is NULL, CELL_FALSE otherwise */
cell_state_t mi0_cmpN (vec128bec_t* v) {
  uint8_t i;

  for(i = 1; i <= 128; i++) {
     if (get_cell(v, i) == CELL_NULL) {
        return (CELL_TRUE);
     }
  }

  return (CELL_FALSE);
}


/* XOR A and B, place result in D     */
void mi0_xor  (cell_proc_t * restrict cp) {
  int i;

  if (verbose_flag) { printf("mi0_xor\n"); }

  if ((mi0_cmpN(cp->A) == CELL_TRUE) || (mi0_cmpN(cp->B) == CELL_TRUE)) { cp_halt(cp, "MI_XOR: NULL Argument(s)");}
  vsetN(cp->D);

  for (i = 1; i <= 128; i++) {
    set_cell(cp->D, i, cxor(get_cell(cp->A, i), get_cell(cp->B,i)));
  }

  if (verbose_flag) { printf("mi0_xor end\n"); }
}


/* LDA - Load vector into register A */
void mi_LDA (cell_proc_t *restrict cp, vec128bec_t* input) {
  vcopy(input, cp->A);
}

/* LDB - Load vector into register B */
void mi_LDB (cell_proc_t *restrict cp, vec128bec_t* input) {
  vcopy(input, cp->B);
}

/* LDC - Load vector into register C */
void mi_LDC (cell_proc_t *restrict cp, vec128bec_t* input) {
  vcopy(input, cp->C);
}

/* LDD - Load vector into register D */
void mi_LDD (cell_proc_t *restrict cp, vec128bec_t* input) {
  vcopy(input, cp->D);
}

/* Rule 30 : x(n+1,i) = x(n,i-1) xor [x(n,i) or x(n,i+1)] */
cell rule_30(cell left, cell middle, cell right){
  return (cxor (left, cor (middle, right)));
}


/* Evaluate source vector with elemental rule, placing result in dest  */
void eval_rule(rule_t r, vec128bec_t* source, vec128bec_t* dest){
  uint8_t col;
  cell left_cell, right_cell, middle_cell;

  for (col = 1; col <= 128; col++){
     left_cell   = (col == 128) ? get_cell(source, 1)    : get_cell(source, col+1);
     right_cell  = (col == 1)   ? get_cell(source, 128)  : get_cell(source, col-1); 
     middle_cell = get_cell(source, col);
            
     set_cell(dest, col, rule_30(left_cell, middle_cell, right_cell) );   
  }  
}

/* Update Time Seed with current time   
   96 bits of clock information are loaded such that long count is in MSB, 
      short count on LSB with and cyclic in center:
   Short Count (32) -> Random(16) -> Cyclic (32) <- Random(16) <- Long Count (32)
   The rest of the bits are seeded with randomness from previous cycle, or 0 on init
 */ 
void mi0_incSDTIME (cell_proc_t *restrict cp){
  int i;
  char* out_str = NULL;
  TimeSeed seed = time_seed();
  
  if (verbose_flag) {  
    printf("MI2_INCSDTIME: %s\n",timeseed_tostr(&seed));
  }   

  if (mi0_cmpN(cp->R) == CELL_TRUE) {
    vsetZ(cp->SDTIME);         // On init, pre-Seed with 0
  } else {
    vcopy(cp->R, cp->SDTIME);  // pre-Seed with previous output for cycle mode
  }                            
  
  /* MSB 1-32 */
  for (i = 1; i <= 32; i++){
    set_cell(cp->SDTIME, i, binary_to_cell ((seed.long_count  >> (i-1)) & 0x01));  
  }

  /* Center 48-80 */
  for (i=1; i<=32; i++){
     set_cell(cp->SDTIME, (48+i), binary_to_cell ((seed.cyclic >> (i-1)) & 0x01));  
  }

  /* LSB 96-128 */
  for (i = 1; i <= 32; i++){
    set_cell(cp->SDTIME, (96+i), binary_to_cell ((seed.short_count >> (i-1)) & 0x01));  

  }  

  if (verbose_flag) {  
    out_str = fmt_vecbe(cp->SDTIME, FMT_VEC_BINARY);
    printf("SDTIME: %s\n",out_str);
    free(out_str);
  }  

}

/* INCPSI - SHA30 of SDTIME, place result in PSI */
void mi2_incPSI (cell_proc_t *restrict cp){

  if (verbose_flag) {  
    printf("MI2_INCPSI\n");
  }   

  vsetN(cp->PSI);
  if (cp->counter_mode == 1) {
     vcopy(cp->SDR30, cp->A);
  } else {
     vcopy(cp->SDTIME, cp->A);
  }
 
  mi2_sha30(cp);
  vmov(cp->D, cp->PSI);
}

/* SHA30 - Use REGA as seed for R30, generate 2 blocks, move second block to REGD
 */
void mi2_sha30 (cell_proc_t *restrict cp){

  if (verbose_flag) {  
    printf("MI2_SHA30\n");
  }   
 /* pushSD(cp); */
  vsetN(cp->D);
  vmov(cp->A, cp->SDR30);
  mi1_incR30(cp);
  mi1_incR30(cp);
  if (cp->counter_mode) {
    vcopy(cp->SDR30, cp->D);
  } else {
    vmov(cp->SDR30, cp->D);
  }
  
 /* popSD(cp);  // Restore Seeds */
}

/* Generate non-deterministic rand, PSI must already exist 
   Take SHA30 of PSI, then XOR that to PSI. Move result to R
 */

void mi2_genR (cell_proc_t *restrict cp) {
  if (verbose_flag) {  
    printf("MI2_GENR\n");
  }   
  vsetN (cp->R);
  vcopy(cp->PSI, cp->A);    // PSI -> A
  mi2_sha30(cp);            // CELEST(PSI) -> D
  vcopy(cp->PSI, cp->A);    // PSI -> A
  vcopy(cp->D, cp->B);      // A:PSI B: CELEST(PSI)
  mi0_xor(cp);              // D: RAND
  vmov(cp->D, cp->R);       // R: PSI XOR CELEST(PSI)
}


/* Machine instruction R30
   Runs R30 for 128 clock cycles to fill R30 with next segment
   RegA : Seed, copied from SDR30, will set center column to 1 if R30 is zero
   RegB : Buffer
   RegD : Output 
   Copies last state to SDR30
 */
void mi1_incR30 (cell_proc_t *restrict cp){
  uint32_t gen;
  uint32_t center_col = 0;
  char* out_str = NULL;

  if (verbose_flag) {  
    printf("MI1_INCR30\n");
  }   

  if (cp_getstate(cp) != CP_IDLE)       { cp_halt (cp, "INCR30: Celproc not idle");  }

  if (mi0_cmpN(cp->SDR30) == CELL_TRUE) { cp_halt (cp, "INCR30: Null Input on SDR30"); }

  vmov(cp->SDR30, cp->A);    // Copy R30 Seed to input 

  vsetN(cp->B);              // Clear Buffer

  cp->NR = cp->A;            // Set row pointers
  cp->CR = cp->B;

  cp_setstate(cp, CP_RUN);

  center_col = (128 / 2);

  if (mi0_cmpZ(cp->A) == CELL_TRUE) {
    set_cell(cp->A, center_col, CELL_TRUE);
  }

  for (gen = 1; gen <= 128; gen++) {
    if (gen & 0x01) {    // Toggle Row pointers every other gen
      cp->NR = cp->B;
      cp->CR = cp->A;
    } else {
      cp->NR = cp->A;
      cp->CR = cp->B;
    }

    if (verbose_flag) {
      out_str = fmt_vecbe(cp->CR, FMT_VEC_BINARY_TEXT);
      printf("[%03d]     %s %c\n", gen, out_str, cell_to_char(get_cell(cp->CR, center_col)));
      free(out_str);
    }
    
    eval_rule(30, cp->CR, cp->NR);
    
    set_cell(cp->D, gen, get_cell(cp->NR, center_col)); //  Copy center cell to D[gen]
  } 

  vmov(cp->CR, cp->SDR30);    /* Update Seed        */
  vmov(cp->D,  cp->R30);      /* Move result to R30 */

  cp_setstate(cp, CP_IDLE);  
}

/* During Time Quantum :
   I/O input transfers occur (A,B loaded), get_neighbor calls resolved 
   TimeSeed (SDTIME) is updated
   R30 is atomically advanced by 128 bits
   PSI is atomically updated
   R is atomically updated based on PSI
   1 Process Cycle occurs
   I/O output transfers occur (D loaded)
   Can be reduced by 1 quantum when incR30 is moved to previous cycle
 */
void mi5_time_quantum(cell_proc_t* restrict cp){
  if (verbose_flag) { printf ("Time Quantum\n");}

  if (!(cp_getstate(cp) == CP_IDLE)) { cp_halt(cp, "time_quantum : CellProc not Idle"); }
  /* Increment time seed if not in counter mode */
  if (cp->counter_mode == 1) {
    if (verbose_flag) { printf ("Time Quantum - counter mode");}
  /* printf("mi5_time_quantum SDR30:\n");
    print_vec(cp->SDR30); */
  } else {
    if (verbose_flag) { printf ("Time Quantum - incSDTIME\n");}
    mi0_incSDTIME(cp);
  }
  
  if (verbose_flag) { printf ("Time Quantum - incR30\n");}
  mi1_incR30(cp);  /* SDR30 -> SDR30 */
  vcopy(cp->SDR30, cp->R);
  if (verbose_flag) { printf ("Time Quantum - incPSI\n");}
  /* mi2_incPSI(cp); */
  if (verbose_flag) { printf ("Time Quantum - genR\n");}
 /* mi2_genR(cp); */
  if (verbose_flag) { printf ("Time Quantum - end\n");}
}

/* Push Seed Registers onto stack */
void pushSD(cell_proc_t* restrict cp) {
  frame_push(cp->SDTIME, cp->Stack);
  frame_push(cp->SDR30,  cp->Stack);
  if (verbose_flag) { 
     printf("\nPUSHSD\n%s\n", frame_to_str(cp->Stack, FMT_VEC_BINARY_TEXT)); 
  }  
}

/* Pop Seed Registers from stack */
void popSD(cell_proc_t* restrict cp){
  if (verbose_flag) { 
     printf("\nPOPSD\n%s\n", frame_to_str(cp->Stack, FMT_VEC_BINARY_TEXT)); 
  }   
  frame_pop(cp->Stack, cp->SDR30);
  frame_pop(cp->Stack, cp->SDTIME);  
}

/* Push General Purpose Registers onto stack */
void pushGP(cell_proc_t* restrict cp) {
  frame_push(cp->A, cp->Stack);
  frame_push(cp->B, cp->Stack);
  frame_push(cp->C, cp->Stack);
  frame_push(cp->D, cp->Stack);

  if (verbose_flag) { 
     printf("PUSHGP\n%s\n", frame_to_str(cp->Stack, FMT_VEC_BINARY_TEXT)); 
  }  
}

/* Pop General Purpose Registers from stack */
void popGP(cell_proc_t* restrict cp) {
  frame_pop(cp->Stack, cp->D);
  frame_pop(cp->Stack, cp->C);
  frame_pop(cp->Stack, cp->B);
  frame_pop(cp->Stack, cp->A);

  if (verbose_flag) { 
     printf("POPGP\n%s\n", frame_to_str(cp->Stack, FMT_VEC_BINARY_TEXT)); 
  }  
}

/* Check clocks
 * Return 0 if OK, print error and halts otherwise
 */
int check_clocks()
{  
  int r = 0;

  if (verbose_flag) { printf("Clock Check\n"); }
   
  /* Test Clock */

  TimeSeed s0 = time_seed();
  TimeSeed s1 = time_seed();

  if (verbose_flag) {
    printf("s0: %s\n",timeseed_tostr(&s0));
    printf("s1: %s\n",timeseed_tostr(&s1));
  }

  /* Error if TimeSeeds did not change */
  if (ts_cmpE(&s0, &s1) == CELL_TRUE) {  
     halt("Clock failed\n");
  }

  return (r);
}




