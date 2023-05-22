#include <fstream>
#include <string.h>
#include <vector>
#include <math.h>
#include <string>
#include <unordered_map>
#include <bitset>
#include <Windows.h>
#include <algorithm>
using namespace std;

/*
* Tan Jing Jie 1804560 P2
* Jacynth Tham Ming Quan 1801600 P2
*/

// for input use
struct location {
	string loc_name;
	int x;
	int y;
};

vector<location> all;
vector<location> city;
vector<location> parcel;
location depot;
int ALL_SIZE, CITY_SIZE, PARCEL_SIZE;

vector<int> city_num;
vector<int> parcel_num;

void get_input() {

	//input
	ifstream input;
	input.open("./tsp.txt");

	location temp;
	while (input >> temp.loc_name) {
		input >> temp.x;
		input >> temp.y;
		switch (temp.loc_name.at(1)) {
		case 'a':
			parcel.push_back(temp);
			break;
		case 'o':
			city.push_back(temp);
			break;
		default:
			depot = temp;
		}
	}
	input.close();

	all.push_back(depot);
	all.insert(all.end(), city.begin(), city.end());
	all.insert(all.end(), parcel.begin(), parcel.end());

	ALL_SIZE = all.size();
	CITY_SIZE = city.size();
	PARCEL_SIZE = parcel.size();
}

void do_output(vector<string> answer) {
	ofstream output("./solution.txt");
	for (string a : answer) {
		output << a << endl;
	}
	output.close();
}

double get_distance(location a, location b) {
	return sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));
}

// held karp ----------------------------------------------

const double INF = numeric_limits<double>::infinity();

// for each possible pair
struct pairset {
	vector<int> city_set; //{1,2,3}
	vector<int> parcel_set; //{4,5,6}
	int total_size;

	pairset(vector<int> c, vector<int> p) {
		city_set = c;
		parcel_set = p;
		total_size = c.size() + p.size();
	}

	pairset() {
		city_set = vector<int>();
		parcel_set = vector<int>();
		total_size = 0;
	}

};


vector<vector<double>> roadmap;

void build_roadmap() {

	for (int i = 0; i < ALL_SIZE; i++) {
		vector<double> row;
		for (int j = 0; j < ALL_SIZE; j++) {
			if (((j > 0 && i > 0 && i <= CITY_SIZE && j > CITY_SIZE)
				|| (i > 0 && j > 0 && i > CITY_SIZE && j <= CITY_SIZE)
				|| (j == 0 && i <= CITY_SIZE)
				|| (i == 0 && j > CITY_SIZE))
				&& (i != j)
				) {
				row.push_back(get_distance(all[i], all[j]));
			}
			else {
				row.push_back(INF);
			}

		}
		roadmap.push_back(row);
	}

	//debug
	/*
	for (int i = 0; i < ALL_SIZE; i++) {
		for (int j = 0; j < ALL_SIZE; j++) {
			cout << roadmap[i][j] << "\t" ;
		}
		cout << endl;
	}
	*/
}

int binomial_coefficient(int n, int k) {
	vector<int> c(k + 1);
	c[0] = 1; // set 1, any nC0 is 1
	for (int i = 1; i <= n; i++){
		for (int j = min(i, k); j > 0; j--)
			c[j] = c[j] + c[j - 1];
	}
	return c[k];
}

int generate_total_calculation(int p, int c) {
	int ans = p * c;

	for (int i = 2; i < c + 1; i++) {
		for (int j = i - 1; j < i + 1; j++) {
			ans += binomial_coefficient(p, i) * binomial_coefficient(c, j);
		}
	}

	return ans;
}

//produce [subset list] =[{},{}]
vector<vector<int>> generate_combinations(vector<int> to_twist, int output_size) { // n C r = to_twist.size C output_size
	vector<vector<int>> combination_list;
	int to_twist_size = to_twist.size();

	//normally no happen
	if (to_twist_size < output_size) {
		return combination_list;
	}

	//initialise
	vector<int> count(output_size);
	vector<int> generated(output_size);
	for (int i = 0; i < output_size; i++) {
		count[i] = i;
		generated[i] = to_twist[i];
	}

	// construt combination
	int difference = to_twist_size - output_size;
	while (true) {
		//pushback (1st loop is push original first few node according to outputsize)
		combination_list.push_back(generated);

		int i = output_size - 1;
		while (i >= 0 && count[i] >= difference + i) {
			i--;
		}

		if (i < 0) {
			break;
		}

		int j = count[i] + 1;
		while (i < output_size) {
			count[i] = j;
			generated[i] = to_twist[j]; //to be push back;
			j++; i++;
		}
	}

	return combination_list;
}

vector<string> compute_path_held_karp() {

	for (int i = 1; i <= CITY_SIZE; i++) {
		city_num.push_back(i);
	}
	for (int i = CITY_SIZE + 1; i < ALL_SIZE; i++) {
		parcel_num.push_back(i);
	}

	// Step 1
	// prepare a dictionary map to save distance for passed route
	// key = all_visited_node__last_visited_node(where u come from)
	// double = commulative distance, long = last_node (where u currently at)
	unordered_map<string, pair<double, unsigned long>> umap;

	int set_size = generate_total_calculation(PARCEL_SIZE, CITY_SIZE);
	vector<pairset> pairset_list(set_size);
	int pairset_index = 0;

	vector<int> routes, combset, nextset;
	string key;
	unsigned long bits, prev;
	bool flag;

	// Symbol convention: each type of bracket is a set {set}  [set list]  (set list list)
	//                                        [pairset(1 set city + 1 set parcel)] (pairset list)

	// Step 2
	// generate a entry in umap from depot to parcel (one to one)
	// generate a pairset from first parcel to city
	// handle case when only have 1 parset - 1 city only
	// input: depot, city {1,2,3}, parcel {4,5,6}
	// output: umap: [1_0,Pair(distance,1),...]
	//         pairset: ([{1},{4}], [{1},{5}], [{1},{6}],[{2},{4}],...,[{3},{6}])

	for (int i : parcel_num) {
		key = to_string(1 << i) + '_' + to_string(i);
		umap[key] = make_pair(roadmap[0][i], 0); //[0][i] is distance of depot to parcel, 0 mean from depot

		vector<int> parcel_set = { i };

		for (int j : city_num) {
			vector<int> city_set = { j };
			pairset_list[pairset_index++] = pairset(city_set, parcel_set);

			// loop and store the route from each parcel to city when only 1 city
			// j is always same in this case
			// detail explaination at following portion of code - generate all possible route
			if (CITY_SIZE == 1) {
				for (int bit : city_set) { // translate to bits indication
					bits |= 1 << bit;
				}
				for (int bit : parcel_set) {
					bits |= 1 << bit;
				}
				routes.push_back(bits); //enter possible routes, to be use at the end (returning depot)
			}
		}
	}



	// Step 3
	
	// generate posible city combination list when 1 city, 2 cities, ..., or n cities are picked 
	// **parcel list is created directly in loop
	// input: city_num = 4
	// output: each row indicate a ilteration of loop
	// ( [{1},{2},{3},{4}]
	//   [{1,2},{1,3},{1,4},{2,3},{2,4},{2,5}...{3,4}]
	//   [{1,2,3},{1,2,4},...,{2,3,4}]
	//   [{1,2,3,4}] )
	vector<vector<vector<int>>> city_set_list_list(CITY_SIZE);
	for (int c = 1; c < CITY_SIZE + 1; c++) {
		city_set_list_list[c - 1] = generate_combinations(city_num, c);
	}

	// to generate all possible route from start parcel to end city
	// input: city set and parcel set
	// output: generate pairset and store in pairset list
	//         possible complete route 
	// route here no have order characteristic, just indicate skiped which parcel or not
	int cz1 = CITY_SIZE + 1; // extract for speed up
	for (int c = 2; c < cz1; c++) {

		// generate possible parcel combination list, c=1,2,3...
		vector<vector<int>> parcel_set_list = generate_combinations(parcel_num, c); //c=2,{1,2,3,4} --> {1,2}, {1,3}, ...,{3,4}

		for (vector<int> parcel_set : parcel_set_list) {
			int k1 = c + 1;// extract for speed up
			for (int k = c - 1; k < k1; k++) {

				//retrieve from the city combination list
				vector<vector<int>> city_set_list = city_set_list_list[k - 1];

				for (vector<int> city_set : city_set_list) {
					pairset_list[pairset_index++] = pairset(city_set, parcel_set);

					//if complete route (mean passed all city and the parcel have same amount of city)
					if ((k + c) == 2 * CITY_SIZE) {
						// {4,5,6}--> 0|(1<<bit) --> 0|10000 --> 10000|100000 -->...-->111 0000 (112)
						// eg: add at 8th in 11011, hence 10011011, add at 7th in 10011011, hence 11011011
						bits = 0;

						for (int bit : city_set) {
							bits |= 1 << bit;
						}
						for (int bit : parcel_set) {
							bits |= 1 << bit;
						}

						routes.push_back(bits); //enter possible routes (to be use for returning depot)
					}
				}
			}
		}
	}



	bits = 0;
	int parent = 0;
	double min = INF;
	double temp_min = 0;
	// Step 4
	// get the distance for all possible complete route
	// input: pairset
	// process: fill the fundamental route from pairset (until completed route as saved in routes) in umap
	//          obtain parent route from umap and add for cummulative
	// output: count shortest commulative distance and store in umap

	// iterate all pair
	// purposely remove a node(visited) and find the min path from umap
	// let 110011 is the route, hence, find 10011 from umap
	// hence the key look like {10011}_?, we are finding the ? indicate the optimum parent(optimum route)
	for (pairset ps : pairset_list) {
		bits = 0;
		combset = vector<int>(ps.total_size);
		int index = 0;

		// get back the all visited node
		// city + concatenate with upper 
		// combine set, so the combined set looks {cccppp}
		for (int bit : ps.city_set) {
			bits |= 1 << bit; // translate to bits indication
			combset[index++] = bit; //put in the existed bit
		}
		for (int bit : ps.parcel_set) {
			bits |= 1 << bit; // translate to bits indication
			combset[index++] = bit; //put in the existed bit
		}

		// prepare for determine which bit to remove
		// eg. city{1,2,3},parcel{4,5} --> city, so nextset=city
		if (ps.city_set.size() == ps.parcel_set.size()) {
			// i need go parcel
			nextset = ps.parcel_set;
		}
		else {
			// i need go city
			nextset = ps.city_set;
		}

		// find the lowest cost for next possible adding node
		// {1,2,3},{4,5,6}
		for (int cs : combset) { //set of combination of city set and parcel set

			prev = bits & ~(1 << cs); //remove a node //cs=3 --> 1011 & ~(1<<3)-->1011& ~(1000)-->1011 &0111 -->>11 -->3 //9-->3,1011->11 

			parent = -1; //start with no parent
			min = INF; // set to max
			flag = false; 

			// find the min path from umap among the key that contain same visited node.
			// the key look like {10011}_?, we are finding the ? indicate the optimum parent(optimum route)
			for (int ns : nextset) {
				// if key is existed
				// find route(except the removed node) from umap
				key = to_string(prev) + '_' + to_string(ns); // previous key //

				if (umap.find(key) != umap.end()) { //cannot find mean no valid route
					flag = true;

					temp_min = umap.at(key).first + roadmap[ns][cs];

					// exchange if find min value
					if (temp_min <= min) {
						min = temp_min;
						parent = ns; //get parent for next layer
					}
				}
			}

			if (flag) {
				key = to_string(bits) + '_' + to_string(cs); //create new entry
				umap[key] = make_pair(min, parent);
			}
		}
	}

	// reinitialize value
	parent = 0;
	bits = 0;
	min = INF;
	temp_min = 0;

	// Step 5
	// loops all possible path from city to depot
	// route bits is identical to umap bits
	// after loop is complete the minimum distance is recorded by parent = k, and route is recorded by bits
	for (int route_bits : routes) {
		for (int k : city_num) {
			key = to_string(route_bits) + '_' + to_string(k); // previous key
			temp_min = umap.at(key).first + roadmap[k][0]; //k = from city to depot

			// exchange value if find more minimum(optimum) distance
			if (temp_min <= min) {
				min = temp_min;
				parent = k;
				bits = route_bits;
			}
		}
	}

	// Step 6
	// the parent variable is now the bottom node
	// {0<--route(candidate_city<--p3<--...<--c1<--p1<--0)}
	// back propagate to find full path
	vector<string> answer(2 * CITY_SIZE + 2, "Depot");

	for (int i = 2 * CITY_SIZE; i > 0; i--) {

		// map path index-element to col list index-element
		answer[i] = all[parent].loc_name; //here parent is from above function set, hence, is actually current node

		key = to_string(bits) + '_' + to_string(parent);
		bits = bits & ~(1 << parent); //remove most left bit
		parent = umap.at(key).second; //update parent for next iteration

	}

	return answer;
}

// pure brute force -------------------------------------------

// for pure brute force method
vector<vector<int>> permute(int n, int r) {
	vector<vector<int>> result;

	if (n == 0)
		return result;

	vector<int> v;
	for (int i = 0; i < n; i++) {
		v.push_back(i);
	}

	do {
		vector<int> temp;
		for (int i = 0; i < r; i++) {
			temp.push_back(v[i]);
		}
		if (!(find(result.begin(), result.end(), temp) != result.end())) {
			result.push_back(temp);
		}
	} while (std::next_permutation(v.begin(), v.end()));

	return result;
}

vector<string> compute_path_pure_brute_force() {

	vector<vector<int>> parcel_list;
	vector<vector<int>> city_list;
	parcel_list = permute(PARCEL_SIZE, CITY_SIZE);
	city_list = permute(CITY_SIZE, CITY_SIZE);

	int PARCEL_LIST_SIZE = parcel_list.size();
	int CITY_LIST_SIZE = city_list.size();
	vector<double> from_begin;
	vector<double> to_end;
	vector<vector<double>> path;

	vector<int> parcel_ans;
	vector<int> city_ans;
	double min_distance = DBL_MAX;

	// calculate distance from depot to each possible parcel
	for (int i = 0; i < PARCEL_SIZE; i++) {
		from_begin.push_back(get_distance(depot, parcel[i]));
	}

	// calculate distance from point to depot
	for (int i = 0; i < CITY_SIZE; i++) {
		to_end.push_back(get_distance(city[i], depot));
	}

	// calculate distance from each path
	for (int i = 0; i < CITY_SIZE; i++) {
		vector<double> path_hor;
		for (int j = 0; j < PARCEL_SIZE; j++) {
			path_hor.push_back(get_distance(city[i], parcel[j]));
		}
		path.push_back(path_hor);
	}

	// loop combination
	for (int i = 0; i < CITY_LIST_SIZE; i++) {
		vector<int> city_c = city_list[i];
		for (int j = 0; j < PARCEL_LIST_SIZE; j++) {
			vector<int> parcel_c = parcel_list[j];

			double distance = from_begin[parcel_c[0]];
			for (int k = 0; k < CITY_SIZE; k++) {
				distance += path[city_c[k]][parcel_c[k]];
				if (k + 1 != CITY_SIZE) {
					distance += path[city_c[k]][parcel_c[k + 1]];
				}
				else {
					distance += to_end[city_c[k]];
				}
			}

			//find minimum distance
			if (distance < min_distance) {
				min_distance = distance;
				parcel_ans = parcel_c;
				city_ans = city_c;
			}

		}
	}

	//prepare answer
	vector<string> answer;

	answer.push_back(depot.loc_name);
	for (int i = 0; i < CITY_SIZE; i++) {
		answer.push_back(parcel[parcel_ans[i]].loc_name);
		answer.push_back(city[city_ans[i]].loc_name);
	}
	answer.push_back(depot.loc_name);

	return answer;
}

// main--------------------------------------------------//

int main() {
	// 0. set process priority
	HANDLE hprocess = GetCurrentProcess();
	SetPriorityClass(hprocess, REALTIME_PRIORITY_CLASS);

	// 1. get input
	get_input();


	// 2. compute process
	vector<string> answer;

	// when held_karp is better than pure brute force
	// however normally will go to held karp due to 7 point only can accept up to 3 red and 3 green
	if (ALL_SIZE > 7) {
		// 2ai. preprocess - create roadmap
		build_roadmap();
		// 2aii. compute path (main algo)
		answer = compute_path_held_karp(); // (2^n)*(n^2)
	}
	else {
		// 2b. compute path (main algo)
		answer = compute_path_pure_brute_force(); //n!
	}

	// 3. output
	do_output(answer);

	return 0;
}