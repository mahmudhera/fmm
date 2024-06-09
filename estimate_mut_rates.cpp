/*
In this script, we will estimate the mutation rates.
Command line arguments: two fasta files, one for a ref genome, and one output file, where
the rates will be written.
*/


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <set>
#include <cassert>

using namespace std;

// Function to read the fasta file
map<string, string> read_fasta(string filename) {
    ifstream file(filename);
    map<string, string> sequences;
    string line;
    string name;
    string seq;
    while (getline(file, line)) {
        if (line[0] == '>') {
            if (!seq.empty()) {
                sequences[name] = seq;
            }
            name = line.substr(1);
            seq = "";
        } else {
            seq += line;
        }
    }
    sequences[name] = seq;
    return sequences;
}


// Function to read the first sequence (by calling the read_fasta function)
string read_first_seq(string filename) {
    map<string, string> sequences = read_fasta(filename);
    return sequences.begin()->second;
}


// Function to get all kmers of a given length from a sequence
vector<string> get_kmers(string seq, int k) {
    vector<string> kmers;
    for (int i = 0; i < seq.size() - k + 1; i++) {
        kmers.push_back(seq.substr(i, k));
    }

    // turn the vector into a set to remove duplicates
    sort(kmers.begin(), kmers.end());
    kmers.erase(unique(kmers.begin(), kmers.end()), kmers.end());

    return kmers;
}


// Function that estimates number of kmers with single substitution and deletion
// returns both estimates
pair<int, int> estimate_single_sub_del(vector<string> kmers_orig, vector<string> kplusone_mers_orig, vector<string> kmers_mut) {
    int num_kmer_single_subst = 0;
    int num_kmer_single_del = 0;

    char bases[] = {'A', 'C', 'G', 'T'};

    // create a hash set of kmers_orig, and kmers_mut, use c++11
    set<string> kmers_orig_set(kmers_orig.begin(), kmers_orig.end());
    set<string> kmers_mut_set(kmers_mut.begin(), kmers_mut.end());

    for (string kmer : kmers_orig) {
        
        // generate all kmers that are 1 subst away from the original kmer
        for (int i = 0; i < kmer.size(); i++) {
            for (char base : bases) {
                
                if (base == kmer[i]) {
                    continue;
                }

                string kmer_subst = kmer;
                kmer_subst[i] = base;

                // if kmer_subst in kmers_orig_set, then continue
                if (kmers_orig_set.find(kmer_subst) != kmers_orig_set.end()) {
                    continue;
                }

                // if kmer_subst in kmers_mut_set, then increment num_kmer_single_subst
                if (kmers_mut_set.find(kmer_subst) != kmers_mut_set.end()) {
                    num_kmer_single_subst++;
                }

            }
        }

    }


    // iterate over all k+1 mers in the original genome
    for (string kplusone_mer : kplusone_mers_orig) {
        
        // generate all kmers that are 1 del away from the original k+1 mer
        for (int i = 0; i < kplusone_mer.size(); i++) {
            string kmer_del = kplusone_mer.substr(0, i) + kplusone_mer.substr(i + 1);

            // if kmer_del in kmers_orig_set, then continue
            if (kmers_orig_set.find(kmer_del) != kmers_orig_set.end()) {
                continue;
            }

            // if kmer_del in kmers_mut_set, then increment num_kmer_single_del
            if (kmers_mut_set.find(kmer_del) != kmers_mut_set.end()) {
                num_kmer_single_del++;
            }

        }

    }

    return make_pair(num_kmer_single_subst, num_kmer_single_del);  
        
}


// function to estimate the mutation rates, returns three doubles
// substitution rate, insertion rate, and deletion rate
tuple<double, double, double> estimate_mut_rates(int len_orig, int len_mut, int num_kmer_single_subst, int num_kmer_single_del, int num_A_orig, int num_A_mut) {
    
    double L = len_orig;
    double L2 = len_mut;
    double fA = num_A_orig;
    double fA_mut = num_A_mut;
    double S = num_kmer_single_subst;
    double D = num_kmer_single_del;

    // print all these values
    cout << "L: " << L << endl;
    cout << "L2: " << L2 << endl;
    cout << "fA: " << fA << endl;
    cout << "fA_mut: " << fA_mut << endl;
    cout << "S: " << S << endl;
    cout << "D: " << D << endl;

    double val1, val2;

    val1 = 3.0 * (fA_mut - 1.0*L2/4.0) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  );
    val2 = 3.0 * (4.0 * S * fA_mut / (4.0*S + 3.0*D) - L2 * S / (4.0*S + 3.0*D) ) / (L - 4.0 * fA);

    cout << val1 << " " << val2 << endl;

    assert(val1 == val2);

    val1 = 3.0 * (- fA + 1.0*L/4) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  );
    val2 = 3.0 * S / (4.0*S + 3.0*D);

    cout << val1 << " " << val2 << endl;

    assert(val1 == val2);

    val1 = 3.0 * (fA_mut - fA + 1.0*L/4 - 1.0*L2/4) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  );
    val2 = 3.0 * (4.0 * S * fA_mut / (4.0*S + 3.0*D) - L2 * S / (4.0*S + 3.0*D) ) / (L - 4.0 * fA) 
                        + 3.0 * S / (4.0*S + 3.0*D);

    cout << val1 << " " << val2 << endl;
    
    // assert that the two values are equal up to 1e-6
    assert(abs(val1 - val2) < 1e-6);
    

    double subst_rate = 3.0 * (fA_mut - fA + 1.0*L/4 - 1.0*L2/4) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  );
    //double subst_rate = 3.0 * (4.0 * S * fA_mut / (4.0*S + 3.0*D) - L2 * S / (4.0*S + 3.0*D) ) / (L - 4.0 * fA) 
    //                    + 3.0 * S / (4.0*S + 3.0*D);
    
    double del_rate = 1.0 * D * subst_rate / S;

    double ins_rate = 1.0 * L2 / L - 1.0 + del_rate;

    return make_tuple(subst_rate, ins_rate, del_rate);
}


// main function
int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " <genome1> <genome2> <output_file>" << endl;
        return 1;
    }

    string genome1_filename = argv[1];
    string genome2_filename = argv[2];
    string output_filename = argv[3];

    string str_orig = read_first_seq(genome1_filename);
    string str_mut = read_first_seq(genome2_filename);

    // get all kmers of length 21
    vector<string> kmers_orig = get_kmers(str_orig, 21);
    vector<string> kmers_mut = get_kmers(str_mut, 21);

    // get all kmers of length 22 in the original genome
    vector<string> kplusone_mers_orig = get_kmers(str_orig, 22);

    // estimate the number of kmers with single substitution and deletion
    pair<int, int> estimates = estimate_single_sub_del(kmers_orig, kplusone_mers_orig, kmers_mut);

    // extract them to variables
    int num_kmer_single_subst = estimates.first;
    int num_kmer_single_del = estimates.second;

    // calculate the mutation rates
    tuple<double, double, double> rates = estimate_mut_rates(str_orig.size(), str_mut.size(), num_kmer_single_subst, num_kmer_single_del, count(str_orig.begin(), str_orig.end(), 'A'), count(str_mut.begin(), str_mut.end(), 'A'));

    // extract the rates
    double subst_rate = get<0>(rates);
    double ins_rate = get<1>(rates);
    double del_rate = get<2>(rates);

    // print the rates
    cout << "Substitution rate: " << subst_rate << endl;
    cout << "Insertion rate: " << ins_rate << endl;
    cout << "Deletion rate: " << del_rate << endl;

    // read the first line of genome2
    string genome2_first_line;
    ifstream file(genome2_filename);
    getline(file, genome2_first_line);

    // format of this line: "> mutated_13110_11715_12488_60639"
    // extract the numbers from this line
    stringstream ss(genome2_first_line);
    string token;
    vector<int> numbers;
    bool first = true;
    while (getline(ss, token, '_')) {
        if (first) {
            first = false;
            continue;
        }
        numbers.push_back(stoi(token));
    }

    cout << "Using the true values...." << endl;

    // print these numbers
    for (int num : numbers) {
        cout << num << " ";
    }

    cout << endl;

    // S, D, I, N_sh = these numbers
    int S = numbers[0];
    int D = numbers[1];
    int I = numbers[2];
    int N_sh = numbers[3];

    // calculate the rate using these
    rates = estimate_mut_rates(str_orig.size(), str_mut.size(), S, D, count(str_orig.begin(), str_orig.end(), 'A'), count(str_mut.begin(), str_mut.end(), 'A'));

    // extract the rates
    subst_rate = get<0>(rates);
    ins_rate = get<1>(rates);
    del_rate = get<2>(rates);

    // print the rates
    cout << "Substitution rate: " << subst_rate << endl;
    cout << "Insertion rate: " << ins_rate << endl;
    cout << "Deletion rate: " << del_rate << endl;

    return 0;
}