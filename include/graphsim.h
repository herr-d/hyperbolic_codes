#ifndef GRAPHSIM_H
#define GRAPHSIM_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>

#include <loccliff.h>
#include <basic_includes.hpp>
#include <cstring>
#include <cstdlib>
#include <unordered_set>

using namespace std;

/*! All vertices in a graph state are numbered beginning with 0. To specify
auch an index, the type VertexIndex (which is just unsigned long) is always
used */
typedef unsigned long VertexIndex;

/*! A GraphRegister object maintains a list of its vertices (qubits), each
described by an object of the class QubitVertex described here.*/
struct QubitVertex {
   /*!byprod is the vertex operator (VOp) associated with the qubit (the name
   stems from the term 'byproduct operator' used for the similar concept in
   the one-way quantum computer.*/
   LocCliffOp byprod;
   /*! neigbors is the adjacency list for this vertex */
   std::unordered_set<VertexIndex> neighbors;
   /*! Upon construction, a qubit vertex is initialised with the Hadamard
   operation as VOp, and with wmpty neighbor list. This makes it represent
   a |0>. */
   QubitVertex (void)
     : byprod (lco_H) {};
};


/*! As we often iterate over sublists of GraphRegister::vertices, this
iterator typedef is a handy abbreviation. */
typedef vector<QubitVertex>::iterator VertexIter;
/*! Another iterator, this one for the adjacency lists QubitVertex::neigbors,
and subsets. */
typedef std::unordered_set<VertexIndex>::iterator VtxIdxIter;

/*! A constant version of VertexIter */
typedef vector<QubitVertex>::const_iterator VertexIterConst;
/*! A constant version of VtxIdxIter */
typedef std::unordered_set<VertexIndex>::const_iterator VtxIdxIterConst;



struct ConnectionInfo {
   bool wasEdge;
   bool non1;
   bool non2;
};

//! A quantum register.
/*! GraphRegister is the central class of graphsim. It represents a register of qubits
that can be entangled with each other. It offers functions to initialize the register,
let gates operate on the qubits, do measurements and print out the state. */
class GraphRegister {
  public:
   /*! This vector stores all the qubits, represented as QubitVertex objects. The index
   of the vector is usually taken as of type VertexIndex. */
   vector<QubitVertex> vertices;
   GraphRegister (VertexIndex numQubits, int randomize = -1);
   GraphRegister (GraphRegister& gr);
   ~GraphRegister () {};
   void local_op (VertexIndex v, LocCliffOp o);
   void hadamard (VertexIndex v);
   void phaserot (VertexIndex v);
   void bitflip (VertexIndex v);
   void inline init(const std::vector<size_type>& v){
      for(auto & i : v)
         vertices[i].byprod = lco_H;  
   }

   void inline clear(){
      for(auto& i : vertices){
         i.byprod = lco_H;
      }
   }

   void phaseflip (VertexIndex v);
   void cphase (VertexIndex v1, VertexIndex v2);
   void cnot (VertexIndex vc, VertexIndex vt);
   int measure (VertexIndex v, LocCliffOp basis = lco_Z,
      bool* determined = NULL, int force = -1);
   void invert_neighborhood (VertexIndex v);
   void print_adj_list (ostream& os = cout) const;
   void print_adj_list_line (ostream& os, VertexIndex i) const;
  private:
   void add_edge (VertexIndex v1, VertexIndex v2);
   void del_edge (VertexIndex v1, VertexIndex v2);
   void toggle_edge (VertexIndex v1, VertexIndex v2);
   int graph_Z_measure (VertexIndex v, int force = -1);
   int graph_Y_measure (VertexIndex v, int force = -1);
   int graph_X_measure (VertexIndex v, bool* determined = NULL, int force = -1);
   void toggle_edges (const std::unordered_set<VertexIndex> vs1,
      const std::unordered_set<VertexIndex> vs2);
   bool remove_byprod_op (VertexIndex v, VertexIndex use_not);
   void cphase_with_table (VertexIndex v1, VertexIndex v2);
   ConnectionInfo getConnectionInfo (VertexIndex v1, VertexIndex v2);
};

#endif
