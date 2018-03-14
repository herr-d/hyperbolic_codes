void printout_code(){
	std::cout << "\n######## X stabilizer ########" << std::endl;
	for(auto stabil : code.X_stabilizer){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## Z stabilizer ########" << std::endl;
	for(auto stabil : code.Z_stabilizer){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## X boundaries ########" << std::endl;
	for(auto stabil : code.X_boundary){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## Z boundaries ########" << std::endl;
	for(auto stabil : code.Z_boundary){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## X operators ########" << std::endl;
	for(auto stabil : code.X_operator){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n######## Z operators ########" << std::endl;
	for(auto stabil : code.Z_operator){
		for(auto qubit : stabil){
			std::cout << qubit << " ";
		}
		std::cout << std::endl;
	}
}