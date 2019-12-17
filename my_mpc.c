#include "my_mpc_helper.c"

void initiate_state(struct State *s){
    s->x= 0; s->y= 4; s->v=2; s->psi= deg2rad(20);
    s->cte= poly_eval(coeffs, s->x) - s->y;
    double psi_des= atan(coeffs[1]+ (2* coeffs[2]* s->x) );
    s->epsi= s->psi- psi_des;
}

int main(){

    double x_arr[20], y_arr[20];
    struct State ss; 
    initiate_state(&ss);
    
    for(int i=0; i<20; i++){
        x_arr[i]= ss.x;
        y_arr[i]= ss.y;
        ss= optimzier(ss);
    }
    
    
    clean_and_print(x_arr, y_arr);
}