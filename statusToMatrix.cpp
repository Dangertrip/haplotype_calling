#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/trim_all.hpp>

using namespace std;



typedef struct haplotype 
{
    int id;
    // string chr;
    // int start;
    int end;
    string mcode; // Methylation code
    unsigned long long mcode_number;

};

vector<map<string, map<int, vector<haplotype>>>> sample_haplo;
vector<map<string, vector<int>>> sample_haplo_cpg_location;


void haplotypeFormat(const std::string & inputFileName)
{
    ifstream haploFile;
    haploFile.open(inputFileName);
    string line;
    int id = 0;
    map<string, map<int, vector<haplotype>>> haplotype_map;
    map<string, vector<int>> cpg_position;
    while (!haploFile.eof())
    {
        //cout << id << "\t";
        getline(haploFile, line);
        boost::trim(line);
        if (line.length()==0) continue;
        vector<string> line_content;
        vector<string> cpg_locations;

        /*
         * Line format: chr1    10000,10002,10004   101
         * which represents: chromsomeName, CpGLocations, MethyaltionStatus
         * Split each line by \t;
         * methylation_code stores the MethylationStatus. "1" is added to the front of the methylation_code(modifiedMethylationCode generated) 
         * for the computational conveniency. 
         */
        boost::split(line_content, line, boost::is_any_of("\t"));
        string chr = line_content[0];
        boost::split(cpg_locations, line_content[1], boost::is_any_of(","));
        string methy_code = line_content[2];
        //for (int i = 0; i<cpg_locations.size();i++) cout << cpg_locations[i] << "\n";
        int end = stoi(cpg_locations.back());

        //cout << "Line load.\t";

        /*
         * For each haplotype, a map of <int startPosition, int modifiedMethylationCode> is generated to capture the overlap parts of two haplotypes. 
         * For example: 
         * If a line is: chr1   10000,10002,10004   101
         * We generate a map contains three haplotypes with the same id:
         * chr1     10000, 10002, 10004     101
         * chr1     10002, 10004    01
         * chr1     10004   1
         * For each sample, haplotype will be splited into multiple parts(N=length) for the next XOR computation.
         */
        if (!haplotype_map.count(chr)) 
        {
            map<int, vector<haplotype>> newchr_map;
            vector<int> cpg_pos;
            haplotype_map[chr] = newchr_map;
            cpg_position[chr] = cpg_pos;

        }
        //cout << "split start.\t";
        for (int i = 0; i < cpg_locations.size()-1;i++) // Overlap region should have at least 2 CpGs.
        {
            //cout << "split" << i << "\t";
            haplotype line_to_haplo;
            //cout << cpg_locations[i] << "\t";
            int start = stoi(cpg_locations[i]);
            line_to_haplo.id = id;
            line_to_haplo.end = end;
            line_to_haplo.mcode = methy_code.substr(i);
            reverse(line_to_haplo.mcode.begin(), line_to_haplo.mcode.end());
            string methy_code_modified = "1" + line_to_haplo.mcode;
            //cout << methy_code_modified << "\n";
            line_to_haplo.mcode_number = stoull(methy_code_modified, nullptr, 2);
            if (haplotype_map[chr].count(start))
            {
                haplotype_map[chr][start].push_back(line_to_haplo);
            }
            else
            {
                cpg_position[chr].push_back(start);
                vector<haplotype> haplo_from_this_cpg;
                haplo_from_this_cpg.push_back(line_to_haplo);
                haplotype_map[chr][start] = haplo_from_this_cpg;
            }
        }
        id++;
        //cout << "split end.\n";
    }
    cout << "Loading finished.\n";
    // Sort all starting position.
    for (map<string, vector<int>>::iterator iter = cpg_position.begin(); iter != cpg_position.end(); iter++)
    {
        sort(iter->second.begin(), iter->second.end());
    }
    cout << "Sorting finished.\n";
    sample_haplo.push_back(haplotype_map);
    sample_haplo_cpg_location.push_back(cpg_position);

};


unsigned long long cantor(int k1, int k2)
{
    unsigned long long a = k1+k2, b = k1+k2+1;
    return a*b/2 + (unsigned long long)k2;
}


int popcount64d(unsigned long long x) // From wikipedia
{
    int count;
    for (count=0; x; count++)
        x &= x - 1;
    return count;
}


float hamming_distance(unsigned long long k1, int k1_length, unsigned long long k2, int k2_length)
{
    int min_length = min(k1_length, k2_length);
    unsigned long long mask = (1 << (min_length - 1)) - 1;
    unsigned long long result = (k1^k2)&mask;
    return (float)popcount64d(result)/(float)min_length;
}



void distanceCompute(vector<string> input_file_names, string output_file_name)
{
    
    int sample_size = input_file_names.size();

    for (int i = 0; i < sample_size; i++)
    {
        cout << "Loading file: " << input_file_names[i] << "\n";
        haplotypeFormat(input_file_names[i]);
    }
    
    float distance[sample_size][sample_size];
    float distance_sum[sample_size][sample_size];
    int distance_count[sample_size][sample_size];

    for (int i = 0; i < sample_size; i++)
    {
        map<string, map<int, vector<haplotype>>>& haplotype_map_a = sample_haplo[i];
        map<string, vector<int>>& cpg_position_a = sample_haplo_cpg_location[i];

        for (int j = i+1; j < sample_size; j++)
        {
            map<string, map<int, vector<haplotype>>>& haplotype_map_b = sample_haplo[j];
            map<string, vector<int>>& cpg_position_b = sample_haplo_cpg_location[j];
            set<unsigned long long> used_pair;
            for (map<string, vector<int>>::iterator iter = cpg_position_a.begin(); iter != cpg_position_a.end(); iter++)
            {
                string chrom = iter->first;
                vector<int>& starting_position_list = iter->second;
                if (!cpg_position_b.count(chrom)) continue;
                for (int start_index = 0; start_index < starting_position_list.size(); start_index++)
                {
                    int cpg_start = starting_position_list[start_index];
                    if (!haplotype_map_b[chrom].count(cpg_start)) continue;
                    vector<haplotype>& haplo_list_a = haplotype_map_a[chrom][cpg_start], haplo_list_b = haplotype_map_b[chrom][cpg_start];
                    for (int ind_a = 0; ind_a<haplo_list_a.size(); ind_a++)
                    {
                        haplotype& haplo_a = haplo_list_a[ind_a];
                        for (int ind_b = 0; ind_b<haplo_list_b.size(); ind_b++)
                        {
                            haplotype& haplo_b = haplo_list_b[ind_b];
                            unsigned long long marker = cantor(haplo_a.id, haplo_b.id);
                            if (used_pair.count(marker)) continue;
                            used_pair.insert(marker);
                            float dist = hamming_distance(haplo_a.mcode_number, haplo_a.mcode.length(), haplo_b.mcode_number, haplo_b.mcode.length());
                            distance_sum[i][j] += dist;
                            distance_count[i][j] += 1;

                        }
                    }
                }
            }

        cout << i << "\t" << j << "\t" << distance_sum[i][j] << "\t" << distance_count[i][j] << "\n";
        }
    }


    vector<string> output_content;
    string header = ","+input_file_names[0];
    for (int i = 1; i < sample_size; i++)
    {
        header = header + "," + input_file_names[i];
    }
    output_content.push_back(header+"\n");
    cout << header;
    for (int i = 0; i < sample_size; i++)
    {
        string output_line = input_file_names[i];
        for (int j = 0; j < sample_size; j++)
        {
            string marker = "NA";
            if (j>i) marker = to_string(distance_sum[i][j]/(float)distance_count[i][j]);
            output_line = output_line + "," + marker;
        }
        cout << output_line;
        output_content.push_back(output_line+"\n");
    }
    ofstream output_stream;
    output_stream.open(output_file_name);
    for (int i = 0; i < output_content.size(); i++)
    {
        cout << output_content[i];
        output_stream << output_content[i];
    }
    output_stream<<flush;
    output_stream.close();

    // Testing prints

    /*cout << sample_haplo.size() << "\n";
    for (int i = 0; i < sample_haplo.size(); i++) {
        for (map<string, map<int, vector<haplotype>>>::iterator iter = sample_haplo[i].begin(); iter != sample_haplo[i].end(); iter++)
        {
            cout << "haplotype count: " << iter->first << "\t" << iter->second.size() << "\n";
        }
    }
    
    cout << sample_haplo_cpg_location.size() << "\n";
    for (int i = 0; i < sample_haplo_cpg_location.size(); i++)
    {
        for (map<string, vector<int>>::iterator iter = sample_haplo_cpg_location[i].begin(); iter != sample_haplo_cpg_location[i].end(); iter++)
        {
            cout << "Start position count: " << iter->first << "\t" << iter->second.size() << "\n";
        }
    }
    
    

    for (int i = 0; i < sample_size; i++)
    {
        for (int j = i+1; j < sample_size; j++)
        {
            map<string, map<int, vector<haplotype>>> haplotype_map_a, haplotype_map_b;
            map<string, vector<int>> cpg_position_a, cpg_position_b;
            
        }
    }*/

};


int main(int argc, char* argv[]) 
{
    string type;
    string db_path;
    string in_file;
    string out_file;
    //string sep("\t");
    int col = -1;

    if (argc > 1) 
    {
        type = string(argv[1]);
    }

    int opt;
    opterr = 0;
    string input_file = "";
    string output_file = "";
    vector<string> haplotype_file_list;

    while ((opt = getopt(argc, argv, "i:o:")) != -1) 
    {
        switch(opt)
        {
            case 'i':
                //printf("option=i, opt=%d, optarg=%s, optind=%d, optopt=%d\n",  opt, optarg, optind, optopt);
                input_file = optarg;
                break;
            case 'o':
                //printf("option=o, opt=%d, optarg=%s, optind=%d, optopt=%d\n",  opt, optarg, optind, optopt);
                output_file = optarg;
                break;
            default:
                break;

        }
    }
    // printf(input_file.c_str());
    // printf(output_file.c_str());
    ifstream input_list_file;
    input_list_file.open(input_file);
    string line;
    if (!input_list_file)
    {
        cout << "Can not open file: " << input_file << " !\n";
        return 1;
    }
    cout << "Reading sample haplotype file list from: " << input_file << "\n";
    ifstream testInputFile;
    while (!input_list_file.eof())
    {
        getline(input_list_file, line);    
        boost::trim(line);
        testInputFile.open(line);
        cout << line << "\n" << line.length() << "\n";
        if (line.length()==0) continue;
        if (!testInputFile)
        {
            cout << "Can not open file: " << line << " !\n";
            return 1;
        }
        testInputFile.close();
        haplotype_file_list.push_back(line);
    }
    input_list_file.close();

    distanceCompute(haplotype_file_list, output_file);
    //for (int i = 0; i<haplotype_file_list.size();i++)
    //{
    //    cout << haplotype_file_list[i] << "\n";
    //}

    return 0;
}


