#ifndef GEN_SURFACE_CODE
#define GEN_SURFACE_CODE
/*

	This file is only used to check for errors in the autotune wrapper library.
	It adds a very simple Surface code generator.

*/

size_type get_data_qubit_id(int x, int y){
	//calculate id
	int id = y/2 * (2*code.r-1);
	id += y%2 * code.r;
	id += x/2;
	return id;
}



void generate_surface_code(){
	code.num_qubits = 2*code.r*code.r - 2*code.r + 1;

	ParityCheck Xbound_one;
	ParityCheck Xbound_two;
	ParityCheck Zbound_one;
	ParityCheck Zbound_two;


	ParityCheck Xdata_one;
	ParityCheck Xdata_two;
	ParityCheck Zdata_one;
	ParityCheck Zdata_two;

	//add stabilizers
	for(int x = 0; x < 2*code.r -1; ++x){
		for(int y = 0; y < 2*code.r -1; ++y){
			if((x%2 == 0 && y%2 == 0) || x%2 == 1 && y%2 ==1){ //data qubit
				if(x == 0){
					code.XBoundaryQubits.insert(get_data_qubit_id(x,y));
				}
				if(x == 2*code.r-2){
					code.XBoundaryQubits.insert(get_data_qubit_id(x,y));
				}
				if(y == 0){
					code.ZBoundaryQubits.insert(get_data_qubit_id(x,y));
				}
				if(y == 2*code.r-2){
					code.ZBoundaryQubits.insert(get_data_qubit_id(x,y));
				}
			}
			else if(x%2 == 1 && y%2 == 0){ // X stabilizer
				ParityCheck tmp;
	
				for(int dx=-1; dx <=1; dx=dx+2)
					if(x+dx >= 0 && x+dx < 2*code.r - 1)
						tmp.insert(get_data_qubit_id(x+dx,y));

					for(int dy=-1; dy <=1; dy=dy+2)
						if(y+dy >= 0 && y+dy < 2*code.r - 1)
							tmp.insert(get_data_qubit_id(x,y+dy));

				code.X_stabilizer.push_back(tmp);
				if(x == 1){
					Xbound_one.insert(code.X_stabilizer.size()-1);
				}
				if(x == 2*code.r-3){
					Xbound_two.insert(code.X_stabilizer.size()-1);
				}
			}

			else{ // Z stabilizer
				ParityCheck tmp;
				for(int dx=-1; dx <=1; dx=dx+2)
					if(x+dx >= 0 && x+dx < 2*code.r - 1)
							tmp.insert(get_data_qubit_id(x+dx,y));

				for(int dy=-1; dy <=1; dy=dy+2)
					if(y+dy >= 0 && y+dy < 2*code.r - 1)
							tmp.insert(get_data_qubit_id(x,y+dy));

				code.Z_stabilizer.push_back(tmp);
				if(y == 1){
					Zbound_one.insert(code.Z_stabilizer.size()-1);
				}
				if(y == 2*code.r-3){
					Zbound_two.insert(code.Z_stabilizer.size()-1);
				}
			}
		}
	}

	//add test
	code.X_stabilizer.back().insert(code.num_qubits);
	++code.num_qubits;



	code.X_boundary.push_back(Xbound_one);
	code.X_boundary.push_back(Xbound_two);
	code.Z_boundary.push_back(Zbound_one);
	code.Z_boundary.push_back(Zbound_two);

	//add logical operators
	ParityCheck X_logical;
	ParityCheck Z_logical;
	for(int i = 0; i < code.r; ++i){
		X_logical.insert(get_data_qubit_id(code.r-1,2*i));
		Z_logical.insert(get_data_qubit_id(2*i,code.r-1));
	}
	code.X_operator.push_back(X_logical);
	code.Z_operator.push_back(Z_logical);

}


#endif //GEN_SURFACE_CODE
