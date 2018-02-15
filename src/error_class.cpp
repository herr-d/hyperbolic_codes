#include <error_class.hpp>

//
// BEGIN BASE CLASS: NO ERRORS
//

void NoError::hadamard(size_type pos)
{
    state_.hadamard(pos);
    return;
}

void NoError::Pauli_X(size_type pos)
{
    state_.local_op(pos,lco_X);
    return;
}

void NoError::Pauli_Z(size_type pos)
{
    state_.local_op(pos,lco_Z);
    return;
}

void NoError::Pauli_Y(size_type pos)
{
    state_.local_op(pos,lco_Y);
    return;
}

int NoError::measurement(size_type pos)
{
    return state_.measure(pos);
}

void NoError::cnot(size_type pos, size_type target)
{
    state_.cnot(pos,target);
    return;
}


// END BASE CLASS

// BEGIN REAL ERROR MODEL



void ErrorClass::hadamard(size_type pos)
{
    NoError::hadamard(pos);
    //TODO
    return;
}

void ErrorClass::Pauli_X(size_type pos)
{
    NoError::Pauli_X(pos);
    //TODO
    return;
}

void ErrorClass::Pauli_Y(size_type pos)
{
    NoError::Pauli_X(pos);
    //TODO
    return;
}

void ErrorClass::Pauli_Z(size_type pos)
{
    NoError::Pauli_X(pos);
    //TODO
    return;
}

int ErrorClass::measurement(size_type pos)
{
    int result = NoError::measurement(pos);
    //TODO
    return result;
}

void ErrorClass::cnot(size_type pos, size_type target)
{
    NoError::cnot(pos,target);
    //TODO
    return;
}
