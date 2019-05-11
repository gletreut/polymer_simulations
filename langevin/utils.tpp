using namespace std;

template<typename T>
void utils::load_matrix(istream &mystream, vector<vector<T> > &cmap){
	string line;
	stringstream convert;
	int i,j,n;
	T val;
	typename
	vector<vector<T> >::iterator it;

	while (getline(mystream,line)){
		convert.clear();
		convert.str(line);

		if ( (convert >> i) && (convert >> j) && (convert >> val) ){
			n=i+1;
			if (j+1 > n) n=j+1;

			// size issues
			if (cmap.size() < n)
				cmap.resize(n);
			for (it=cmap.begin(); it!=cmap.end(); it++){
				if (it->size() < n)
					it->resize(n);
			}

			// put value in table
			cmap[i][j]=val;
		}
	}

	cout << "MATRIX WITH SIZE " << n << " x " << n << " IMPORTED" << endl;

}

template<typename T>
void utils::print_matrix(ostream &mystream, vector<vector<T> > &cmap){
	mystream << left;
	//mystream << fixed << setprecision(8);
	mystream << scientific << setprecision(8);

	for (size_t i=0; i<cmap.size(); i++){
		for (size_t j=0; j<cmap.size(); j++){
			mystream << setw(10) << i;
			mystream << setw(10) << j;
			mystream << setw(20) << cmap[i][j];
			mystream << endl;
		}
	}
}

template<typename T>
void utils::load_map(istream &mystream, map<string,T> &params){
	string line, key;
	stringstream convert;
	T val;
	cout << left;

	while (getline(mystream, line)){
		convert.clear();
		convert.str(line);

		if (convert >> key) {// there is a string
			if (convert >> val) {// there is a value
				params[key]=val;
				cout << setw(30) << key;
				cout << val << endl;
			}
		}
	}

	return;
}

template<typename T>
void utils::print_map(ostream &mystream, map<string,T> &params){
	int lmax,l;
	typename
	map<string,T>::iterator it;

	lmax=0;
	for (it=params.begin(); it!=params.end(); it++){
		l=it->first.length();
		if (l>lmax) lmax=l;
	}

	lmax+=2;

	mystream << left;
	for (it=params.begin(); it!=params.end(); it++){
		mystream << setw(lmax) << it->first.c_str();
		mystream << it->second;
		mystream << endl;
	}
	return;
}

