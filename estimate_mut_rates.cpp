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


    // placeholder
    num_kmer_single_del = -1;

    return make_pair(num_kmer_single_subst, num_kmer_single_del);  
        
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

    // write the estimates to the output file
    ofstream output_file(output_filename);
    output_file << estimates.first << " " << estimates.second << endl;
    output_file.close();

    return 0;
}