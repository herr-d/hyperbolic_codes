/*
 * This function prints out the code object in a readable format.
 * Intended to use for debugging.
 *
*/
void printout_code(){
	std::cout << "\n######## X stabilizer ######## size: " << code.X_stabilizer.size() << std::endl;
	for(auto stabil : code.X_stabilizer){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## Z stabilizer ######## size: " << code.Z_stabilizer.size() << std::endl;
	for(auto stabil : code.Z_stabilizer){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## X boundaries ######## size: " << code.X_boundary.size()<< std::endl;
	for(auto stabil : code.X_boundary){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## Z boundaries ######## size: " << code.Z_boundary.size() << std::endl;
	for(auto stabil : code.Z_boundary){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## X operators ######## size:" << code.X_operator.size() << std::endl;
	for(auto stabil : code.X_operator){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## Z operators ######## size:" << code.Z_operator.size() << std::endl;
	for(auto stabil : code.Z_operator){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}
}