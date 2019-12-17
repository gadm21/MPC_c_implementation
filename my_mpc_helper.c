#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

//optimizer variables
const int horizon_steps= 3;
double motor_dof[6] = {0, 1, 2, 3, 4, 5};  //6
double steering_dof[9] = { -20, -10, -5, -2, 0, 2, 5, 10, 20};  //9

//model variables
double coeffs[3]= {0, 0 ,1}; //1x2+ 0x+ 0
const int voltage2speed_weight = 2;
const double dt= 0.1;

//cost variables
const double cte_weight = 9000;
const double epsi_weight = 800;
const double v_weight = 50;
const double acceleration_cost_weight = 600;
const double steering_cost_weight = 100;
const double change_steer_cost_weight = 100;
const double change_motor_cost_weight = 100;
const double ref_cte = 0;
const double ref_epsi = 0;
const double ref_v = 1 * voltage2speed_weight;

#define M_PI 3.14159265358979323846
#define ui unsigned int

//helper functions
double deg2rad(double deg){ return deg * M_PI / 180; }
double poly_eval(double coeff[3], double x){
    double result= 0.0;
    result+= (coeff[2]* x*x)+
             (coeff[1]* x)+
             coeff[0];
    return result;
}

struct State{
    double x, y, v, psi, cte, epsi;
    double MV, steering_angle;
};

struct State model(struct State *s, double mv, double steering_angle){

    double new_v= voltage2speed_weight * mv;
    double new_psi= s->psi+ (new_v * deg2rad(steering_angle) * dt);
    double new_x= s->x+ (cos(new_psi)* new_v* dt);
    double new_y= s->y+ (sin(new_psi)* new_v* dt);
    double new_cte= poly_eval(coeffs, new_x) - new_y;
    double psi_des= atan(coeffs[1]+ (2* coeffs[2]* new_x) );
    double new_epsi= new_psi- psi_des;

    struct State new_s;
    new_s.x= new_x; new_s.y= new_y; new_s.epsi= new_epsi;
    new_s.steering_angle= steering_angle; new_s.MV= mv; 
    new_s.psi= new_psi; new_s.cte= new_cte; new_s.v= new_v;
    return new_s;
}

double cost(struct State s0, struct State s1, struct State s2, struct State s3){

    double result=0.0;

    result+= (cte_weight * abs(s0.cte - ref_cte) +
              cte_weight * abs(s1.cte - ref_cte) + 
              cte_weight * abs(s2.cte - ref_cte)+
              cte_weight * abs(s3.cte - ref_cte));

    result+= (epsi_weight * abs(s0.epsi - ref_epsi) +
             epsi_weight * abs(s1.epsi - ref_epsi) +
             epsi_weight * abs(s2.epsi - ref_epsi) +
             epsi_weight * abs(s3.epsi - ref_epsi));          

    
    result+= (v_weight * abs(s0.v - ref_v) +
              v_weight * abs(s1.v - ref_v) +
              v_weight * abs(s2.v - ref_v) +
              v_weight * abs(s3.v - ref_v));
              
    result+= (steering_cost_weight * abs(s0.steering_angle) +
              steering_cost_weight * abs(s1.steering_angle) +
              steering_cost_weight * abs(s2.steering_angle));
    
    //didn't add acceleration_cost_weight because, unlike a real vehicle, in an RC car you need MV to keep the speed constant
    

    result+= (change_steer_cost_weight * abs(s1.steering_angle - s0.steering_angle) +
              change_steer_cost_weight * abs(s3.steering_angle - s2.steering_angle));
    
    result+= (change_motor_cost_weight * abs(s1.MV - s0.MV) +
              change_motor_cost_weight * abs(s3.MV - s2.MV));
              
    return result;
}


struct State optimzier(struct State s){

    int motor_dof_size= 6;//sizeof(motor_dof)/ sizeof(motor_dof[0]);
    int steering_dof_size= 9;//sizeof(steering_dof) / sizeof(steering_dof[0]);

    int step1_size= motor_dof_size* steering_dof_size;
    int step2_size= step1_size* step1_size;
    int step3_size= step2_size* step1_size;
    
    struct State step1_states[step1_size];
    struct State step2_states[step2_size];
    


    for(int i=0; i<motor_dof_size; i++)
        for(int j=0; j<steering_dof_size; j++)
            step1_states[(i*steering_dof_size)+j]= model(&s, motor_dof[i], steering_dof[j]);
    
    for (int step=0; step< step1_size; step++)
        for(int i=0; i<motor_dof_size; i++)
            for(int j=0; j<steering_dof_size; j++)  
                step2_states[(i*steering_dof_size)+j]= model(&step1_states[step], motor_dof[i], steering_dof[j]);
    
    struct State last_state;
    int lowest_cost_id= 0;
    int lowest_cost= INT_MAX;
    int id=0;
    for (int step=0; step< step2_size; step++)
        for(int i=0; i<motor_dof_size; i++)
            for(int j=0; j<steering_dof_size; j++){
                
                struct State final_state= model(& step2_states[step], motor_dof[i], steering_dof[j]);

                double current_cost= cost(s,step1_states[step/step1_size], step2_states[step], final_state);
                if(current_cost< lowest_cost){
                    lowest_cost_id= id;//(step* step1_size)+ (i* steering_dof_size)+j;
                    lowest_cost= current_cost;
                }                

                id++;
            }
    
    struct State next_state= step1_states[lowest_cost_id/ step2_size];
    return next_state;

}



void clean_and_print(double arrx[20], double arry[20]){
    printf("plt.plot([ ");
    for(int i=0; i<20; i++) if(i==19) printf("],[ \n"); else printf("%.2f, ", arrx[i]);
    for(int i=0; i<20; i++) if(i==19) printf("], 'r')\n"); else printf("%.2f, ", arry[i]);
}


void report(struct State s){
    printf("state report-> v:%lf, x:%lf, y:%lf, psi:%lf, vm:%lf, steering:%lf, epsi:%lf, cte:%lf \n",
            s.v, s.x, s.y, s.psi, s.MV, s.steering_angle, s.epsi, s.cte);
}

void full_report(struct State * path){
    int size= sizeof(path)/ sizeof(path[0]);
    for(int i=0; i<size; i++) report(path[i]);
}

