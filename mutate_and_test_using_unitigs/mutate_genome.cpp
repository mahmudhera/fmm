#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <unordered_map>
#include <tuple>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <set>


using namespace std;

std::vector<char> alphabet = {'A', 'C', 'G', 'T'};

std::unordered_map<char, std::vector<char>> substitute_dic = {
    {'A', {'C', 'G', 'T'}},
    {'C', {'A', 'G', 'T'}},
    {'G', {'C', 'A', 'T'}},
    {'T', {'C', 'G', 'A'}}
};

class mutation_model {
public:
    mutation_model(int seed, size_t orig_length, double p_s, double p_d, double d)
        : seed(seed), orig_length(orig_length), p_s(p_s), p_d(p_d), d(d) {
        srand(seed);
    }

    // constructor with given original string
    mutation_model(int seed, string orig_string, double p_s, double p_d, double d)
        : seed(seed), orig_string(orig_string), p_s(p_s), p_d(p_d), d(d) {
        srand(seed);
        orig_length = orig_string.length();
    }

    string generate_random_string(const vector<double>& frequency = {}) {
        size_t len = orig_length;
        string orig_string;

        orig_string.reserve(len);
        if (frequency.empty()) {
            for (size_t i = 0; i < len; i++) {
                orig_string.push_back(alphabet[std::rand() % alphabet.size()]);
            }
        } else {
            while( orig_string.length() < len ) {
                double r = static_cast<double>(rand()) / RAND_MAX;
                double cumulative_prob = 0.0;

                for (int j = 0; j < alphabet.size(); j++) {
                    cumulative_prob += frequency[j];
                    if (r < cumulative_prob) {
                        orig_string.push_back(alphabet[j]);
                        break;
                    }
                }
            }
        }

        orig_string.shrink_to_fit();
        this->orig_string = orig_string;

        assert(orig_string.length() == orig_length);

        return orig_string;
    }

    tuple<string, size_t, size_t, size_t, size_t> mutate_string(int k = -1) {
        string new_string = string(this->orig_string);
        vector<int> actions_chosen(orig_length);
        vector<int> insert_lengths(orig_length);
        vector<string> strings_to_insert(orig_length);

        // Choose actions for each position in the string
        for (size_t i = 0; i < orig_length; i++) {
            actions_chosen[i] = choose_action();
            insert_lengths[i] = choose_insert_length();
            strings_to_insert[i] = generate_random_string(insert_lengths[i]);
        }

        // Apply chosen actions to mutate the string
        for (size_t i = 0; i < orig_length; i++) {
            if (actions_chosen[i] == 2) {
                char orig_char = new_string[i];
                new_string[i] = substitute_dic[orig_char][std::rand() % 3];
            } else if (actions_chosen[i] == 3) {
                new_string[i] = ' ';
            }
        }

        // Combine the mutated string and inserted substrings
        std::string final_str = combine_strings(new_string, strings_to_insert);

        // Calculate mutation counts if k is provided
        if (k == -1) {
            return std::make_tuple(final_str, 0, 0, 0, 0);
        }

        size_t num_kmers_single_insertion = 0;
        size_t num_kmers_single_deletion = 0;
        size_t num_kmers_single_substitution = 0;
        size_t num_kmers_no_mutation = 0;

        // sets for kmers
        set<pair<string, string>> kmers_with_single_insertion;
        set<pair<string, string>> kmers_with_single_deletion;
        set<pair<string, string>> kmers_with_single_substitution;

        for (size_t i = 0; i < orig_length - k + 1; i++) {
            int sum1 = 0;
            int sum2 = 0;
            for (int j = i; j < i + k; j++) {
                sum1 += actions_chosen[j];
                sum2 += insert_lengths[j];
            }
            if (sum1 == 3 && sum2 == 0) {
                num_kmers_single_deletion++;
                string corresponding_mut_kmer = "";
                for (int ii = 0; ii < k; ii++) {
                    if (new_string[ii] != ' ') {
                        corresponding_mut_kmer = corresponding_mut_kmer + new_string[ii];
                    }
                    corresponding_mut_kmer += corresponding_mut_kmer + strings_to_insert[ii];
                }
                kmers_with_single_deletion.insert(make_pair(this->orig_string.substr(i, k), corresponding_mut_kmer));
            }
            if (sum1 == 2 && sum2 == 0) {
                num_kmers_single_substitution++;
                string corresponding_mut_kmer = "";
                for (int ii = 0; ii < k; ii++) {
                    if (new_string[ii] != ' ') {
                        corresponding_mut_kmer = corresponding_mut_kmer + new_string[ii];
                    }
                    corresponding_mut_kmer += corresponding_mut_kmer + strings_to_insert[ii];
                }
                kmers_with_single_substitution.insert(make_pair(this->orig_string.substr(i, k), corresponding_mut_kmer));
            }
            if (sum1 == 0 && sum2 == 1) {
                num_kmers_single_insertion++;
                string corresponding_mut_kmer = "";
                for (int ii = 0; ii < k; ii++) {
                    if (new_string[ii] != ' ') {
                        corresponding_mut_kmer = corresponding_mut_kmer + new_string[ii];
                    }
                    corresponding_mut_kmer += corresponding_mut_kmer + strings_to_insert[ii];
                }
                kmers_with_single_insertion.insert(make_pair(this->orig_string.substr(i, k), corresponding_mut_kmer));
            }
            if (sum1 == 0 && sum2 == 0) {
                num_kmers_no_mutation++;
            }
        }


        // open three files and write the kmers to them
        ofstream kmer_file_insertion("kmers_single_insertion");
        ofstream kmer_file_deletion("kmers_single_deletion");
        ofstream kmer_file_substitution("kmers_single_substitution");

        for (auto kmer : kmers_with_single_insertion) {
            kmer_file_insertion << kmer << endl;
        }

        for (auto kmer : kmers_with_single_deletion) {
            kmer_file_deletion << kmer << endl;
        }

        for (auto kmer : kmers_with_single_substitution) {
            kmer_file_substitution << kmer << endl;
        }

        kmer_file_insertion.close();
        kmer_file_deletion.close();
        kmer_file_substitution.close();

        return make_tuple(final_str, num_kmers_single_substitution, num_kmers_single_insertion, num_kmers_single_deletion, num_kmers_no_mutation);
    }

    void test() {
        cout << this->orig_string << endl;
    }

private:
    int seed;
    size_t orig_length;
    double p_s;
    double p_d;
    double d;
    string orig_string;

    int choose_action() {
        double rand_num = static_cast<double>(rand()) / RAND_MAX;
        if (rand_num < p_s) {
            return 2; // Substitution
        } else if (rand_num < p_s + p_d) {
            return 3; // Deletion
        } else {
            return 0; // Stay (no mutation)
        }
    }

    // Choose an insertion length following a geometric distribution
    int choose_insert_length() {
        double rand_num = static_cast<double>(rand()) / RAND_MAX;
        int len = static_cast<int>(log(rand_num) / log(d/(1.0+d)) );
        if (len < 0) len = 0;
        return len;
    }

    // Generate a random string of a given length
    string generate_random_string(int length) {
        string random_str;
        random_str.reserve(length);
        for (int i = 0; i < length; i++) {
            random_str.push_back(alphabet[rand() % alphabet.size()]);
        }
        return random_str;
    }

    // Combine the mutated string with inserted substrings
    string combine_strings(const string& base, const vector<string>& inserts) {
        string final_str;
        final_str.reserve(base.size() + inserts.size());
        for (int i = 0; i < orig_length; i++) {
            if (base[i] != ' ')
                final_str += base[i];
            final_str += inserts[i];
        }
        return final_str;
    }
};

int main(int argc, char* argv[])
{
    // command line arguments: seed, original_filename, p_s, p_d, d, output_filename, kmer_size

    // print usage
    if (argc != 8) {
        cout << "Usage: " << argv[0] << " seed orig_genome_filename p_s p_d d output_filename kmer_size" << endl;
        return 1;
    }

    // parse command line arguments
    int seed = atoi(argv[1]);
    string orig_filename = argv[2];
    double p_s = atof(argv[3]);
    double p_d = atof(argv[4]);
    double d = atof(argv[5]);
    string output_filename = argv[6];
    int ksize = atoi(argv[7]);

    // read original genome fasta file, skip the first line
    ifstream orig_file(orig_filename);
    string line;
    string orig_string = "";
    getline(orig_file, line);
    while (getline(orig_file, line)) {
        orig_string += line;
    }

    // make the orig string all uppercase
    transform(orig_string.begin(), orig_string.end(), orig_string.begin(), ::toupper);

    // strip the original string of all non-ACGT characters
    orig_string.erase(remove_if(orig_string.begin(), orig_string.end(), [](char c) { return c != 'A' && c != 'C' && c != 'G' && c != 'T'; }), orig_string.end());
    
    // create mutation model
    mutation_model model(seed, orig_string, p_s, p_d, d);

    // generate mutated string by calling mutate_string()
    auto mut_res = model.mutate_string(ksize);
    size_t S = get<1>(mut_res);
    size_t I = get<2>(mut_res);
    size_t D = get<3>(mut_res);
    size_t N_shared = get<4>(mut_res);

    // generate output header with following format: "> mutated_S_D_I_N_shared"
    string output_header = "> mutated_" + to_string(S) + "_" + to_string(D) + "_" + to_string(I) + "_" + to_string(N_shared);

    // write header string, then mutated string to file
    ofstream output_file(output_filename);
    output_file << output_header << endl;
    output_file << get<0>(mut_res) << endl;
    output_file.close();

    return 0;
}
