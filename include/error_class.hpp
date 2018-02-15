#ifndef ERRORCLASS_HPP
#define ERRORCLASS_HPP


#include "graphsim.h"
#include <random>


class NoError
{
protected:
    std::mt19937_64 rng;
    GraphRegister& state_;
public:
    NoError(uint64_t seed, GraphRegister& state) : rng(seed), state_(state){
    }

    void Pauli_X(size_type pos);
    void Pauli_Y(size_type pos);
    void Pauli_Z(size_type pos);
    void cnot(size_type pos, size_type target);
    void hadamard(size_type pos);
    void init(size_type pos);
    int measurement(size_type pos);
};


class ErrorClass : NoError{
    ErrorClass(uint64_t seed, GraphRegister& state) : NoError(seed,state){}
    void Pauli_X(size_type pos);
    void Pauli_Y(size_type pos);
    void Pauli_Z(size_type pos);
    void cnot(size_type pos, size_type target);
    void hadamard(size_type pos);
    void init(size_type pos);
    int measurement(size_type pos);
};


#endif
