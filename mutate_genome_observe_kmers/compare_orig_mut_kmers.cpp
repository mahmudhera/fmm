# include <iostream>
# include <vector>
# include <string>
# include <random>
# include <unordered_map>
# include <tuple>
# include <fstream>
# include <algorithm>
# include <cassert>
# include <unordered_set>


using namespace std;

std::vector<char> alphabet = {'A', 'C', 'G', 'T'};

// given a string and kmer size, return its kmers
vector<string> get_kmers(string& s, int ksize) {
    vector<string> kmers;
    for (int i = 0; i < s.length() - ksize + 1; i++) {
        kmers.push_back(s.substr(i, ksize));
    }
    return kmers;
}

// given two vectors of kmers, return the number of shared kmers
size_t get_num_shared_kmers(vector<string>& kmers1, vector<string>& kmers2) {
    unordered_set<string> kmer_set1(kmers1.begin(), kmers1.end());
    unordered_set<string> kmer_set2(kmers2.begin(), kmers2.end());
    int num_shared_kmers = 0;
    for (auto kmer : kmer_set1) {
        if (kmer_set2.find(kmer) != kmer_set2.end()) {
            num_shared_kmers++;
        }
    }
    return num_shared_kmers;
}

// given two vectors of kmers, return the number of kmers with single subst, delt, insert, and shared
tuple<size_t, size_t, size_t, size_t> get_num_kmers_with_single_subst_delt_insert_shared(vector<string>& kmers_orig, vector<string>& kmers_mutated) {
    // create set of kmers from vectors
    unordered_set<string> kmer_set_orig(kmers_orig.begin(), kmers_orig.end());
    unordered_set<string> kmer_set_mutated(kmers_mutated.begin(), kmers_mutated.end());

    // create set of shared kmers
    unordered_set<string> shared_kmer_set;
    for (auto kmer : kmer_set_orig) {
        if (kmer_set_mutated.find(kmer) != kmer_set_mutated.end()) {
            shared_kmer_set.insert(kmer);
        }
    }

    // create a set of k-1 mers from kmers_orig
    unordered_set<string> k_minus_1_mers_set_orig;
    for (auto kmer : kmer_set_orig) {
        int k = kmer.length();
        k_minus_1_mers_set_orig.insert( kmer.substr(0, k-1) );
        k_minus_1_mers_set_orig.insert( kmer.substr(1, k-1) );
    }

    // create a set of k-1 mers from kmers_mutated
    unordered_set<string> k_minus_1_mers_set_mutated;
    for (auto kmer : kmer_set_mutated) {
        int k = kmer.length();
        k_minus_1_mers_set_mutated.insert( kmer.substr(0, k-1) );
        k_minus_1_mers_set_mutated.insert( kmer.substr(1, k-1) );
    }

    // create a hash set of all k-1 -mers in the mutated string that is 1 deletion away from a k-mer in the mutated string
    unordered_set<string> k_minus_1_mers_set_mutated_1_deletion;
    for (auto kmer : kmer_set_mutated) {
        int k = kmer.length();
        for (int i = 0; i < k; i++) {
            string k_minus_1_mer = kmer.substr(0, i) + kmer.substr(i+1, k-i-1);
            k_minus_1_mers_set_mutated_1_deletion.insert(k_minus_1_mer);
        }
    }

    // variables to return
    size_t num_kmers_with_single_subst = 0;
    size_t num_kmers_with_single_delt = 0;
    size_t num_kmers_with_single_insert = 0;
    size_t num_shared_kmers = shared_kmer_set.size();

    // iterate through kmers_orig, check if it has a 1-delt away k-1 mer in k_minus_1_mers_set_mutated_1_deletion
    // if yes, and if it is not in shared set, then it is a kmer with single substitution
    for (auto kmer : kmer_set_orig) {
        int k = kmer.length();
        for (int i = 0; i < k; i++) {
            string k_minus_1_mer = kmer.substr(0, i) + kmer.substr(i+1, k-i-1);
            if (k_minus_1_mers_set_mutated_1_deletion.find(k_minus_1_mer) != k_minus_1_mers_set_mutated_1_deletion.end()) {
                if (shared_kmer_set.find(kmer) == shared_kmer_set.end()) {
                    num_kmers_with_single_subst++;
                }
            }
        }
    }

    // iterate through kmers_orig, check if it has a 1-delt away k-1 mer in k_minus_1_mers_set_mutated
    // if yes, and if it is not in shared set, then it is a kmer with single deletion
    for (auto kmer : kmer_set_orig) {
        int k = kmer.length();
        for (int i = 0; i < k; i++) {
            string k_minus_1_mer = kmer.substr(0, i) + kmer.substr(i+1, k-i-1);
            if (k_minus_1_mers_set_mutated.find(k_minus_1_mer) != k_minus_1_mers_set_mutated.end()) {
                if (shared_kmer_set.find(kmer) == shared_kmer_set.end()) {
                    num_kmers_with_single_delt++;
                }
            }
        }
    }

    // iterate through kmers_mutated, check if it has a 1-delt away k-1 mer in k_minus_1_mers_set_orig
    // if yes, and if it is not in shared set, then it is a kmer with single insertion
    for (auto kmer : kmer_set_mutated) {
        int k = kmer.length();
        for (int i = 0; i < k; i++) {
            string k_minus_1_mer = kmer.substr(0, i) + kmer.substr(i+1, k-i-1);
            if (k_minus_1_mers_set_orig.find(k_minus_1_mer) != k_minus_1_mers_set_orig.end()) {
                if (shared_kmer_set.find(kmer) == shared_kmer_set.end()) {
                    num_kmers_with_single_insert++;
                }
            }
        }
    }

    return make_tuple(num_kmers_with_single_subst, num_kmers_with_single_delt, num_kmers_with_single_insert, num_shared_kmers);

}

int main(int argc, char* argv[]) {
    // ecoli.fasta: original genome name
    // mutation rates: [0.01, 0.05, 0.1]
    // kmer size: 21
    // number of simulations: 10, seed = 0, 1, ..., 9
    // mutated_genome_name = ecoli.fasta_mutated_<subst_rate>_<delt_rate>_<insert_rate>_<seed>.fasta

    // read original genome
    string orig_genome_name = "ecoli.fasta";
    string orig_genome;
    ifstream orig_file(orig_genome_name);
    string line;
    getline(orig_file, line); // skip the first line
    while (getline(orig_file, line)) {
        orig_genome += line;
    }

    // close file
    orig_file.close();

    // p_s, p_d, d = 0.01, 0.05, 0.1, seed = 0, 1, ..., 9
    vector<double> subst_rates = {0.01, 0.05, 0.1};
    vector<double> delt_rates = {0.01, 0.05, 0.1};
    vector<double> insert_rates = {0.01, 0.05, 0.1};
    int ksize = 21;
    int num_simulations = 10;

    // iterate p_s, p_d, d, seed 
    for (double p_s : subst_rates) {
        for (double p_d : delt_rates) {
            for (double d : insert_rates) {
                for (int seed = 0; seed < num_simulations; seed++) {
                    // generate mutated genome filename
                    string mutated_genome_name = orig_genome_name + "_mutated_" + to_string(p_s) + "_" + to_string(p_d) + "_" + to_string(d) + "_" + to_string(seed) + ".fasta";

                    // read mutated genome string, first line format: "> mutated_S_D_I_N_shared", read S, D, I, N_shared
                    ifstream mutated_file(mutated_genome_name);
                    getline(mutated_file, line);
                    size_t S, D, I, N_shared;
                    sscanf(line.c_str(), "> mutated_%zu_%zu_%zu_%zu", &S, &D, &I, &N_shared);

                    // read mutated genome string
                    string mutated_genome;
                    while (getline(mutated_file, line)) {
                        mutated_genome += line;
                    }

                    // close file
                    mutated_file.close();

                    // generate original kmers
                    vector<string> orig_kmers = get_kmers(orig_genome, ksize);

                    // generate mutated kmers
                    vector<string> mutated_kmers = get_kmers(mutated_genome, ksize);

                    // get number of kmers with single subst, delt, insert, and shared
                    auto result = get_num_kmers_with_single_subst_delt_insert_shared(orig_kmers, mutated_kmers);
                    size_t num_kmers_with_single_subst = get<0>(result);
                    size_t num_kmers_with_single_delt = get<1>(result);
                    size_t num_kmers_with_single_insert = get<2>(result);
                    size_t num_shared_kmers = get<3>(result);

                    // print to stdout
                    cout << p_s << " " << p_d << " " << d << " " << seed << " " << S << " " << D << " " << I << " " << N_shared << " " << num_kmers_with_single_subst << " " << num_kmers_with_single_delt << " " << num_kmers_with_single_insert << " " << num_shared_kmers << endl;

                }
            }
        }
    }

    return 0;
}