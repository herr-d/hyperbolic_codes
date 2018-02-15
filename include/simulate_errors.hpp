#ifndef SIMULATE_ERRORS
#define SIMULATE_ERRORS
#include <error_class.hpp>


class Simulator{
private:
    NoError& state_;

    //one cycle of stabilzier measurements
    void cycle_();

    //performs the X-stabilizer on all qubits simultaneously
    void xstabilizer_();

    //performs the Y-stabilizer on all qubits simultaneously
    void zstabilizer_();

public:
    Simulator(NoError &state) : state_(state){};
    void clear();

    //runs the simulation with the given errors

    void simulate();
};

#endif
//SIMULATE_ERRORS