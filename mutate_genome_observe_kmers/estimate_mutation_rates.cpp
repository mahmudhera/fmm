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
pair<uint64_t, uint64_t> estimate_single_sub_del(vector<string> kmers_orig, vector<string> kplusone_mers_orig, vector<string> kmers_mut) {
    uint64_t num_kmer_single_subst = 0;
    uint64_t num_kmer_single_del = 0;

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

    assert(abs(val1 - val2) < 1e-6);

    val1 = 3.0 * (- fA + 1.0*L/4) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  );
    val2 = 3.0 * S / (4.0*S + 3.0*D);

    cout << val1 << " " << val2 << endl;

    assert(abs(val1 - val2) < 1e-6);

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


tuple<double, double, double> estimate_rates_for_pair_of_files_by_kmers(string filename1, string filename2, int k)
{
    string str_orig = read_first_seq(filename1);
    string str_mut = read_first_seq(filename2);

    // get all kmers
    vector<string> kmers_orig = get_kmers(str_orig, k);
    vector<string> kmers_mut = get_kmers(str_mut, k);

    // get all kmers of length k+1 in the original genome
    vector<string> kplusone_mers_orig = get_kmers(str_orig, k+1);

    // estimate the number of kmers with single substitution and deletion
    pair<uint64_t, uint64_t> estimates = estimate_single_sub_del(kmers_orig, kplusone_mers_orig, kmers_mut);

    // extract them to variables
    uint64_t num_kmer_single_subst = estimates.first;
    uint64_t num_kmer_single_del = estimates.second;

    // calculate the mutation rates
    tuple<double, double, double> rates = estimate_mut_rates(str_orig.size(), str_mut.size(), num_kmer_single_subst, num_kmer_single_del, count(str_orig.begin(), str_orig.end(), 'A'), count(str_mut.begin(), str_mut.end(), 'A'));

    return rates;
}


tuple<double, double, double> estimate_rates_for_pair_of_files_by_known_values(string filename1, string filename2, int k)
{
    string str_orig = read_first_seq(filename1);
    string str_mut = read_first_seq(filename2);

    // read the first line of genome2
    string genome2_first_line;
    ifstream file(filename2);
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

    // S, D, I, N_sh = these numbers
    int S = numbers[0];
    int D = numbers[1];
    int I = numbers[2];
    int N_sh = numbers[3];

    // calculate the rate using these
    tuple<double, double, double> rates = estimate_mut_rates(str_orig.size(), str_mut.size(), S, D, count(str_orig.begin(), str_orig.end(), 'A'), count(str_mut.begin(), str_mut.end(), 'A'));

    return rates;
}


// main function
int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " <genome_name> <num_files> <kmer_len>" << endl;
        return 1;
    }

    string genome_filename = argv[1];
    int num_files = stoi(argv[2]);
    int k = stoi(argv[3]);

    // create a vector of mutation rates, keeping the rates as strings
    vector<string> mutation_rates = {"0.01", "0.05", "0.1"};

    for (int i = 0; i < num_files; i++) {
        for (string subst_rate : mutation_rates) {
            for (string del_rate : mutation_rates) {
                for (string ins_rate : mutation_rates) {
                    // use two decimal places for the rates
                    string filename = genome_filename + "mutated_" + subst_rate + "_" + del_rate + "_" + ins_rate + "_" + to_string(k) + "_" + to_string(i) + ".fasta";
                    cout << filename << endl;
                    tuple<double, double, double> rates_by_known_values = estimate_rates_for_pair_of_files_by_known_values(genome_filename, filename, k);
                    tuple<double, double, double> rates_by_kmers = estimate_rates_for_pair_of_files_by_kmers(genome_filename, filename, k);
                    // print: p_s, p_d, p_i, p_s_est1, p_d_est1, p_i_est1, p_s_est2, p_d_est2, p_i_est2
                    cout << subst_rate << " " << del_rate << " " << ins_rate << " ";
                    cout << get<0>(rates_by_known_values) << " " << get<2>(rates_by_known_values) << " " << get<1>(rates_by_known_values) << " ";
                    cout << get<0>(rates_by_kmers) << " " << get<2>(rates_by_kmers) << " " << get<1>(rates_by_kmers) << endl; 
                }
            }
        }
    }

    return 0;
}