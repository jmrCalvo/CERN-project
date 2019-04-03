#include <array>
#include <vector>
#include <random>
#include <map>
#include <set>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
 #include <ctime> 

#include <range/v3/core.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/algorithm/sort.hpp>

#include <generate.h>

namespace {

template <std::size_t N>
struct get_n {
  template <typename T>
  auto operator()(T&& t) const ->
    decltype(std::get<N>(std::forward<T>(t))) {
    return std::get<N>(std::forward<T>(t));
  }
};

const float PI = 3.1415927;

using namespace ranges;

using std::array;
using std::map;
using std::set;
using std::vector;
using std::string;
using std::fstream;
using std::ctime;

}
array<array<float, 3>, 3> pmt_rotation(float, float)
{
  return {{{1.f, 0.f, 0.f},
           {0.f, 1.f, 0.f},
           {0.f, 0.f, 1.f}}};
}

array<float, 3> rotate(const array<array<float, 3>, 3>& rotation,
                       const array<float, 3>& pmt)
{
  return pmt;
}



std::string toBinary(size_t n)
{
    std::string r;
    while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
		while(r.size()!=32){r="0"+r;}
    return r;
}

int main() {

	int num_samples;
	map<int, array<float, 3>> PMTs; //PMT position
	map<int, array<float, 3>> PMTd; //PMT direction
	map<int, array<int, 3>> PMTmetadata; //PMT direction

	map<std::string,int> Probabilities;
  // Part 1: parse km3net_reference.detx into storage
  string line;
  std::ifstream myfile ("../challenge/km3net_reference.detx");
	if (myfile.is_open()){
		int count_lines=0;
    bool FLAG_NEW_DATA=false;
		double DOM_ID=0;double MOD_ID=0;double num_of_DOM=0;
		while(getline (myfile,line)){
			count_lines++;
			if(count_lines>4){//the three first lines are ignored
				if(count_lines==4){//the information about the optical modules
					num_samples=atoi(line.c_str());
					std::cout<<"line 4"<<std::endl;
				}
				else{
					if(count_lines==5 || FLAG_NEW_DATA){//the header of the data detector
						std::stringstream temp;	
						temp<<line;
						FLAG_NEW_DATA=false;
						temp>>DOM_ID;
						temp>>MOD_ID;
						temp>>num_of_DOM;
					}
					else{
						if(line==""){FLAG_NEW_DATA=true;}//is a blanquet line which means there is new data
						else{ //the data stored by the detector 
							std::stringstream temp;	
							temp<<line;
							int PMT_ID=0;float Xpos,Ypos,Zpos; float Xdir,Ydir,Zdir;
							temp>>PMT_ID;
							temp>>Xpos;temp>>Ypos;temp>>Zpos;
							temp>>Xdir;temp>>Ydir;temp>>Zdir;
							//the calibration is not necessary
							//position
							std::array<float,3> arrayTmpPos={Xpos,Ypos,Zpos};
							PMTs[PMT_ID]=arrayTmpPos;
							//direction
							std::array<float,3> arrayTmpDir={Xdir,Ydir,Zdir};
							PMTd[PMT_ID]=arrayTmpDir;
							//metadata
							std::array<int,3> arrayTmpMetadata={DOM_ID,MOD_ID,num_of_DOM};
							PMTmetadata[PMT_ID]=arrayTmpMetadata;

							std::string id=std::to_string(sqrt(pow(Xpos,2)+pow(Ypos,2)+pow(Zpos,2)));
							if(Probabilities.find(id)!=Probabilities.end()){
							Probabilities[id]++;
							}
							else{
							Probabilities[id]=1;
							}
						}
					}
				}	
			}	
    }
    myfile.close();
  }
std::cout<<std::endl<<"made PART 1"<<std::endl;

//for (std::map<std::string,int>::iterator it=Probabilities.begin(); it!=Probabilities.end(); ++it)
  //  std::cout << it->first << " => " << it->second << '\n';
  // Example storage of PMT postions
  //map<int, array<float, 3>> PMTs;

  // Generate random hits
  // The rates are Hz per number of hits on neighbouring PMTs; this
  // needn't be changed.
  const array<float, 4> background_rates{7000.f, 700.f, 70.f, 7.f};
  Generators generators{41431, 2340, background_rates};

  // KM3NeT data arrives in time slices of 100ms. Here this data is
  // generated to a reasonable approximation of reality. How that
  // works is outside the scope of this challenge. Here we look at the
  // equivalent of 1 second of data.
  const size_t n_slices = 10;

  // Select 10 sets of hits per slice
  const size_t n_candidates = 10;

  // Generate slices of 100 ms (1e8 ns)
  const long dt = std::lround(1e8);

  for (size_t i = 0; i < n_slices; ++i) {
    // Values are three numbers packed into a single 32bit integer
    // - 8 bits: time over threshold of this hit (ToT)
    // - 5 bits: PMT_ID
    // - rest  : 100 * (DOM_ID + 1) + MOD_ID + 1
    // The ToT time and ToT are not relevant for this challenge
    auto [times, values] = generate(i * dt, (i + 1) * dt, generators, true);

    // Sort by time, this ensures there is some spread in PMTs within
    // clusters
    ranges::sort(view::zip(times, values), std::less<>{}, get_n<0>{});
		
				
		//auto it2=times.begin();
	map<int,int> Total_Hits;
	
	for(auto it=values.begin();it<values.end();it++){		
			//std::cout<<*(it2)<<"			";
			//std::cout<<*(it)<<"	";
			//size_t ToT=*(it)>>24;

			//size_t PMT_ID_Arr=(*(it)&16252928)>>19;
			//00000000000001111111111111111111 = 524287 100 == 1100100
			//size_t rest=*(it)&524287;
			//size_t rest_DOM=(rest-1)/100;
			//size_t rest_MOD=(rest-1)%100;
			//std::cout<<ToT<<"	"<<PMT_ID_Arr<<"	"<<rest_DOM<<"	"<<rest_MOD<<std::endl;//esto tiene que ser los primero bits
			//it2++;

			//00000000111111111111111111111111 = 16777215
			size_t restyPMT=*(it)&16777215;	
			
			//std::string clave= std::to_string(rest_DOM)+std::to_string(rest_MOD);
			
			if(Total_Hits.find(restyPMT)!=Total_Hits.end()){
					Total_Hits[restyPMT]++;
			}
			else{	
					Total_Hits[restyPMT]=1;
			}

		}

		int means=20;
		int mu=5;
	
	set<int> candidate_Hits;
		for(auto it=Total_Hits.begin();it!=Total_Hits.end();it++){
			if(it->second <(means+mu/2) && it->second >(means-mu/2) && rand()%2!=0){
					candidate_Hits.insert(it->first);
			}
		}
		std::cout<<"Second Part Made"<<std::endl;

    // Part 2: randomly select hits to use as candidates
    // Suggestion: use a normally distributed number of hits with a mean
    // of 20 an a width of 5, but no less than 5

/*
    for (size_t candidate = 0; candidate < n_candidates; ++candidate) {
      size_t hit_start = 0, hit_end = 0;

      // Part 3: rotate selected hits over a grid of directions
      vector<array<float, 3>> result(hit_end - hit_start);

      float d_tp = 5.f * 2.f * PI / (360.f * 360.f);
      for (float theta = 0.f; theta < PI; theta += d_tp) {
        for (float phi = 0.f; phi < 2 * PI; phi += d_tp) {
          for (size_t i = hit_start; i < hit_end; ++i) {

            // Obtain PMT ID
            size_t pmt_id = 0;

            // Calculate rotation
            auto rotation = pmt_rotation(theta, phi);
            // Apply rotation
            result[i] = rotate(rotation, PMTs[pmt_id]);
          }
        }
      }
    }
*/

	for (auto it=candidate_Hits.begin();it!=candidate_Hits.end();it++){
			size_t PMT_ID_Arr=(*(it)&16252928)>>19;
			//00000000000001111111111111111111 = 524287 100 == 1100100
			size_t rest=*(it)&524287;
			size_t rest_DOM=(rest-1)/100;
			size_t rest_MOD=(rest-1)%100;
			std::cout<<PMT_ID_Arr<<"	"<<rest_DOM<<"	"<<rest_MOD<<std::endl;
			
	}


  }
}
