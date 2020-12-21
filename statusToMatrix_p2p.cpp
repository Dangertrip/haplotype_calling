#include <algorithm>
#include <functional>
#include <cctype>
#include <ctime>
#include <locale>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <unistd.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <hash_map>
#include <unordered_set>
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



//vector<map<string, map<int, vector<haplotype>>>> sample_haplo;
//vector<map<string, vector<int>>> sample_haplo_cpg_location;

bool print_detail = false;
unordered_map<string, unordered_map<int, unordered_map<int, vector<haplotype>>>> haplo_hashmap; // hash_map<chromosome, hash_map<start_position, hash_map<cpg_length, vector<haplotype>>>>


void haplotypeFormat(int id, const std::string & inputFileName) // id: file id for the further accumulation
{
    
    time_t init_time = time(0); 


    FILE* haploFile = fopen(inputFileName.c_str(), "rb");
    //string content;
    fseek(haploFile, 0, SEEK_END);
    long size = ftell(haploFile);
    //content.resize(size);
    char *content_arr = (char*)malloc(size+1);

    fseek(haploFile, 0, SEEK_SET);
    int len = fread(content_arr, 1, size, haploFile);
    fclose(haploFile);
    time_t read_time = time(0);

    // For efficiency, we need to deal with the char array here instead of using boost to split the string.
    char *t;
    int tab_count = 0, comma_count = 0;
    int head = 0, word_len = 0;
    string chr;
    int start;
    int pos = 0;

    for (t = content_arr; *t!='\0'; t++)
    {
        //cout << *t;
        if (*t=='\n')
        {
            if (tab_count==2) 
            {
                string methy_code(&content_arr[head], pos-head);
                //cout << "code:  " << "\n";
                //cout << pos-head << "\n";
                //cout << methy_code << "\n";
                
                haplotype line_to_haplo;
                line_to_haplo.id = id;
                line_to_haplo.end = 0;
                line_to_haplo.mcode = methy_code;
                string methy_code_modified = "1" + line_to_haplo.mcode;
                line_to_haplo.mcode_number = stoull(methy_code_modified, nullptr, 2);

                if (haplo_hashmap.find(chr)==haplo_hashmap.end())
                {
                    unordered_map<int, unordered_map<int, vector<haplotype>>> newchr_map;
                    haplo_hashmap[chr] = newchr_map;
                }
                if (haplo_hashmap[chr].find(start)==haplo_hashmap[chr].end())
                {
                    unordered_map<int, vector<haplotype>> newstart_map;
                    haplo_hashmap[chr][start] = newstart_map;
                }
                int mcode_length = methy_code.length();
                if (haplo_hashmap[chr][start].find(mcode_length)==haplo_hashmap[chr][start].end())
                {
                    vector<haplotype> newlength_vector;
                    haplo_hashmap[chr][start][mcode_length] = newlength_vector;
                }
            
                haplo_hashmap[chr][start][mcode_length].push_back(line_to_haplo);
                //cout << chr << "\t" << start << "\t" << mcode_length << "\t" << haplo_hashmap[chr][start][mcode_length].size() << "\n";
                tab_count = 0;
                comma_count = 0;
                head = pos + 1;
            }
            else break;
            
        }
        else if (*t=='\t')
        {
            if (comma_count>0)
            {
                comma_count = 0;
                head = pos + 1;
            } 
            else 
            {
                if (tab_count == 0)
                {
                    chr = string(&content_arr[head], pos-head); // establish string using [head, pos)
                    //cout << "chrom:  " << "\n";
                    //cout << &content_arr[head] << "\n";
                    //cout << pos-head << "\n";
                    //cout << chrom << "\n";
                    //chr = chrom;
                    head = pos + 1;
                }
                if (tab_count == 1)
                {
                    string st(&content_arr[head], pos-head);
                    start = stoi(st);
                    head = pos + 1;
                }
                //else
                //{
                //    string code(&content_arr[head], pos-head);
                //    cout << "code:  " << "\n";
                //    //cout << &content_arr[head] << "\n";
                //    cout << pos-head << "\n";
                //    cout << code << "\n";
                //    methy_code = &code;
                //    head = pos + 1;
                //}
            }
            tab_count++;
        }
        else if (*t==',')
        {
            if (comma_count==0)
            {
                string st(&content_arr[head], pos-head);
                start = stoi(st);
                //cout << "start\t" << start << endl; 
            }
            comma_count++;
        }
        pos++;
    }

    //delete chr;
    //delete t;
    //delete methy_code;


    /*------------------------------------------Code discarded: Old version of data loading-----------------------------------------------
    string content(content_arr);
    vector<string> sections;
    time_t read_time = time(0);
    boost::split(sections, content, boost::is_any_of("\t\n"));
    
    int i = 0;
    int total_sections = sections.size();
    total_sections = total_sections - total_sections%3;

    time_t boost_time = time(0);

    while (i<total_sections)
    {
        //if (i-3>=0 && sections[i+2]==sections[i-1] && sections[i+1]==sections[i-2] && sections[i]==sections[i-3]) continue;
        string chr = sections[i];
        vector<string> cpg_locations;
        //boost::split(cpg_locations, sections[i+1], boost::is_any_of(","));
        string methy_code = sections[i+2];
        int start = stoi(sections[i+1].substr(0, sections[i+1].find(",")));
        int end = 0;
        //int start = stoi(cpg_locations[0]);
        //int end = stoi(cpg_locations.back());

        if (haplo_hashmap.find(chr)==haplo_hashmap.end())
        {
            unordered_map<int, unordered_map<int, vector<haplotype>>> newchr_map;
            haplo_hashmap[chr] = newchr_map;
        }
        if (haplo_hashmap[chr].find(start)==haplo_hashmap[chr].end())
        {
            unordered_map<int, vector<haplotype>> newstart_map;
            haplo_hashmap[chr][start] = newstart_map;
        }
        int mcode_length = methy_code.length();
        if (haplo_hashmap[chr][start].find(mcode_length)==haplo_hashmap[chr][start].end())
        {
            vector<haplotype> newlength_vector;
            haplo_hashmap[chr][start][mcode_length] = newlength_vector;
        }
        haplotype line_to_haplo;
        line_to_haplo.id = id;
        line_to_haplo.end = end;
        line_to_haplo.mcode = methy_code;
        string methy_code_modified = "1" + line_to_haplo.mcode;
        line_to_haplo.mcode_number = stoull(methy_code_modified, nullptr, 2);

        haplo_hashmap[chr][start][mcode_length].push_back(line_to_haplo);

        i = i+3;
    } */
    time_t process_time = time(0);
    cout << "Time diff(read, boost, process)(s): " << difftime(read_time, init_time) << "\t"  << difftime(process_time, read_time) << endl;
/*    
    // Given reading using getline, which is super slow. So this part of the program is discarded.
    while (!haploFile.eof())
    {
        //cout << id << "\t";
        getline(haploFile, line);
        boost::trim(line);
        if (line.length()==0) continue;
        vector<string> line_content;
        vector<string> cpg_locations;

        //
        // Line format: chr1    10000,10002,10004   101
        // which represents: chromsomeName, CpGLocations, MethyaltionStatus
        // Split each line by \t;
        // methylation_code stores the MethylationStatus. "1" is added to the front of the methylation_code(modifiedMethylationCode generated) 
        // for the computational conveniency. 
        //
        boost::split(line_content, line, boost::is_any_of("\t"));
        string chr = line_content[0];
        boost::split(cpg_locations, line_content[1], boost::is_any_of(","));
        string methy_code = line_content[2];
        //for (int i = 0; i<cpg_locations.size();i++) cout << cpg_locations[i] << "\n";
        int start = stoi(cpg_locations[0]);
        int end = stoi(cpg_locations.back());

        //cout << "Line load.\t";

        //
        // For each haplotype, only store itself instead of all sub-haplotypes. 
        //
        if (!haplo_hashmap.count(chr)) 
        {
            unordered_map<int, unordered_map<int, vector<haplotype>>> newchr_map;
            haplo_hashmap[chr] = newchr_map;
        }
        if (!haplo_hashmap[chr].count(start))
        {
            unordered_map<int, vector<haplotype>> newstart_map;
            haplo_hashmap[chr][start] = newstart_map;
        }
        int mcode_length = methy_code.length();
        if (!haplo_hashmap[chr][start].count(mcode_length))
        {
            vector<haplotype> newlength_vector;
            haplo_hashmap[chr][start][mcode_length] = newlength_vector;
        }
        haplotype line_to_haplo;
        line_to_haplo.id = id;
        line_to_haplo.end = end;
        line_to_haplo.mcode = methy_code;
        string methy_code_modified = "1" + line_to_haplo.mcode;
        line_to_haplo.mcode_number = stoull(methy_code_modified, nullptr, 2);

        haplo_hashmap[chr][start][mcode_length].push_back(line_to_haplo);

        //cout << "split start.\t";
        //cout << "split end.\n";
    }
    */

    cout << "Loading finished.\n";
    // Sort all starting position.
    //for (map<string, vector<int>>::iterator iter = cpg_position.begin(); iter != cpg_position.end(); iter++)
    //{
    //    sort(iter->second.begin(), iter->second.end());
    //}
    //cout << "Sorting finished.\n";
    //sample_haplo.push_back(haplotype_map);
    //sample_haplo_cpg_location.push_back(cpg_position);

}

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
    //int min_length = min(k1_length, k2_length);
    //unsigned long long mask = (1 << (min_length - 1)) - 1;
    unsigned long long result = k1^k2;
    return (float)popcount64d(result)/(float)k1_length;
}



void distanceCompute(vector<string> input_file_names, string output_file_name)
{
    
    int sample_size = input_file_names.size();
    
    time_t program_init = time(0);

    for (int i = 0; i < sample_size; i++)
    {
        cout << "Loading file: " << input_file_names[i] << "\n";
        haplotypeFormat( i, input_file_names[i] );
    }

    time_t loading_now = time(0);
    cout << "Time(loading)(s): " << difftime(loading_now, program_init) << endl;
   
    float distance[sample_size][sample_size];
    float distance_sum[sample_size][sample_size] = {0};
    int distance_count[sample_size][sample_size] = {0};

    unordered_map<string, unordered_map<int, unordered_map<int, vector<haplotype>>>>::iterator chr_iter;
    unordered_map<int, unordered_map<int, vector<haplotype>>>::iterator start_index_iter;
    unordered_map<int, vector<haplotype>>::iterator haplo_length_iter;

    for (int i = 0; i<sample_size; i++)
    {
        for (int j = 0; j<sample_size; j++)
        {
            distance_sum[i][j] = 0;
            distance_count[i][j] = 0;
            //cout << i << "\t" << j << "\t" << distance_sum[i][j] << "\t" << distance_count[i][j] << "\n";
        }
    }

    int region_id = 0;

    for (chr_iter = haplo_hashmap.begin(); chr_iter != haplo_hashmap.end(); chr_iter++)
    {
        string chr = chr_iter->first;
        unordered_map<int, unordered_map<int, vector<haplotype>>> haplo_map_chr = chr_iter->second;
        for (start_index_iter = haplo_map_chr.begin(); start_index_iter != haplo_map_chr.end(); start_index_iter++ )
        {
            int start_index = start_index_iter->first;
            unordered_map<int, vector<haplotype>> haplo_map_start = start_index_iter->second;
            for (haplo_length_iter = haplo_map_start.begin(); haplo_length_iter != haplo_map_start.end(); haplo_length_iter++)
            {
                int haplo_length = haplo_length_iter->first;
                vector<haplotype> haplo_vector = haplo_length_iter->second;
                for (int i = 0; i<haplo_vector.size(); i++)
                {
                    for (int j = i+1; j<haplo_vector.size(); j++)
                    {
                        // compute hamming distance for haplotype_i and haplotype_j;
                        haplotype *haplo_a = &haplo_vector[i];
                        haplotype *haplo_b = &haplo_vector[j];
                        if (haplo_a->id==haplo_b->id) continue;
                        float dist = hamming_distance(haplo_a->mcode_number, haplo_a->mcode.length(), haplo_b->mcode_number, haplo_b->mcode.length());
                        if (print_detail){
                        cout << "Cellpair: " << region_id << "\t" << dist << "\t" << haplo_a->id << "\t" << haplo_a->mcode << "\t" << haplo_b->id << "\t" << haplo_b->mcode << "\t" << chr << "\t" << start_index << "\n";
                        }
			//cout << "p2p:\t" << dist << "\t" << haplo_a->id << "\t" << haplo_a->mcode << "\t" << haplo_b->id << "\t" << haplo_b->mcode << "\t" << chr << "\t" << start_index << "\n";
                        distance_sum[haplo_a->id][haplo_b->id] += dist;
                        distance_count[haplo_a->id][haplo_b->id] += 1;
                        distance_sum[haplo_b->id][haplo_a->id] += dist;
                        distance_count[haplo_b->id][haplo_a->id] += 1;
                    }
                }
                region_id++;
            }
        }
    }


    /*
    //==========================================discarded=============================================================
    for (int i = 0; i < sample_size; i++)
    {
        map<string, map<int, vector<haplotype>>>& haplotype_map_a = sample_haplo[i];
        map<string, vector<int>>& cpg_position_a = sample_haplo_cpg_location[i];

        for (int j = i+1; j < sample_size; j++)
        {
            distance_sum[i][j] = 0;
            distance_count[i][j] = 0;
            map<string, map<int, vector<haplotype>>>& haplotype_map_b = sample_haplo[j];
            map<string, vector<int>>& cpg_position_b = sample_haplo_cpg_location[j];
            // set<unsigned long long> used_pair;
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
                            if (haplo_a.mcode.length()!=haplo_b.mcode.length()) continue;
                            // unsigned long long marker = cantor(haplo_a.id, haplo_b.id);
                            // if (used_pair.count(marker)) continue;
                            // used_pair.insert(marker);
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
    */

    for (int i = 0; i<sample_size; i++)
    {
        for (int j = 0; j<sample_size; j++)
        {
            cout << i << "\t" << j << "\t" << distance_sum[i][j] << "\t" << distance_count[i][j] << "\n";
        }
    }
    time_t compute_now = time(0);
    cout << "Time(computation)(s): " << difftime(compute_now, loading_now) << endl;
 
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
            string marker = "-1";
            if (j!=i && distance_count[i][j]>0) 
            {
                marker = to_string(distance_sum[i][j]/(float)distance_count[i][j]);
            }
            else
            {
                marker = "-1";
            }
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

    while ((opt = getopt(argc, argv, "i:o:d:")) != -1) 
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
            case 'd':
                print_detail = true;
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
        // cout << line << "\n" << line.length() << "\n";
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
