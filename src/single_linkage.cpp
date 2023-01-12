#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <math.h>
#include <cstdlib>

using namespace std;

struct Hit {
  string name;
  string target;
  unsigned long begin;
  unsigned long end;
  string orientation;
  double score;
  Hit() { begin = end = 0; }
  friend bool operator< (const Hit& hit1, const Hit& hit2) {
    if (hit1.begin < hit2.begin) {
      return 1;
    }
    return 0;
  }
  friend bool operator>= (const Hit& hit1, const Hit& hit2) {
    if (hit1.begin >= hit2.begin) {
      return 1;
    }
    return 0;
  }
};

typedef set<string> StringSet;
typedef vector<Hit> HitVector;
typedef map<string, HitVector> StringHitVectorMap;

bool Overlap(long cluster_begin, long cluster_end, Hit& h, double dist);
void PrintHit(Hit hit);

int main(int argc, char* argv[]) {
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " feature_file single_linkage_distance\n";
    exit(1);
  }
  string hf = argv[1];
  double dist = atof(argv[2]);

  // read list of hits
  // group hits by chromosome
  ifstream hitfile(hf.c_str());
  assert(hitfile);
  StringHitVectorMap hits;
  StringSet targets;
  StringSet orientations;
  while (!hitfile.eof()) {
    Hit tmphit;
    hitfile >> tmphit.name;
    if (!hitfile.eof()) {
      hitfile >> tmphit.target >> tmphit.orientation;
      hitfile >> tmphit.begin >> tmphit.end >> tmphit.score;
      targets.insert(tmphit.target);
      orientations.insert(tmphit.orientation);
      string tmptarget = (tmphit.target + tmphit.orientation);
      StringHitVectorMap::iterator pos = hits.find(tmptarget);
      if (pos != hits.end()) {
        pos->second.push_back(tmphit);
      }
      else {
        HitVector newhitvector;
        newhitvector.push_back(tmphit);
        hits.insert(pair<string, HitVector>(tmptarget, newhitvector));
      }
    }
  }
  hitfile.close();
  cerr << "Read feature file" << endl;
  int clusterID = 0;
  for (StringSet::iterator pos = targets.begin(); pos != targets.end(); pos++) {
    cerr << "Working on features from " << *pos << endl;
    // separate strands
    for (StringSet::iterator orPtr = orientations.begin(); orPtr != orientations.end(); orPtr++) {
      // sort hits by begins
      string tmptarget = (*pos) + (*orPtr);
      if (hits.find(tmptarget) != hits.end()) {
        HitVector::iterator begin = hits.find(tmptarget)->second.begin(), end = hits.find(tmptarget)->second.end();
        sort(begin, end);
        // now construct clusters
        long cluster_begin = -1, cluster_end = -1;
        string cluster_orientation;
        HitVector current_cluster;
        HitVector empty_vector;
        string target;
        for (HitVector::iterator ptr = begin; ptr != end; ptr++) {
          if (!Overlap(cluster_begin, cluster_end, *ptr, dist)) {
            // new cluster
            if (current_cluster.size() > 0) {
              target = current_cluster.begin()->target;
              cout << "Cluster: slc" << clusterID << " " << tmptarget << " " <<
                cluster_begin << ".." << cluster_end << endl;
              for (HitVector::iterator tmp = current_cluster.begin(); tmp != current_cluster.end(); tmp++) {
                PrintHit(*tmp);
              }
              cout << endl;
            }
            // update clusterID
            clusterID++;
            cluster_begin = (*ptr).begin;
            cluster_end = (*ptr).end;
            current_cluster = empty_vector;
            current_cluster.push_back(*ptr);
          }
          else {
            // update cluster end
            if (cluster_end < (*ptr).end) {
              cluster_end = (*ptr).end;
            }
            current_cluster.push_back(*ptr);
          }
        }
        // print last cluster
        if (current_cluster.size() > 0) {
          target = current_cluster.begin()->target;
          cout << "Cluster: slc" << clusterID << " " << target << " " <<
            cluster_begin << ".." << cluster_end << endl;
          for (HitVector::iterator tmp = current_cluster.begin(); tmp != current_cluster.end(); tmp++) {
            PrintHit(*tmp);
          }
          cout << endl;
        }
      }
    }
  }
}

bool Overlap(long cluster_begin, long cluster_end, Hit& h, double dist) {
  if (cluster_begin >= 0 && h.begin >= 0) {
    if ((cluster_begin <= h.end && h.begin <= cluster_end) ||
      fabs(double(cluster_begin) - double(h.begin)) <= dist ||
      fabs(double(cluster_begin) - double(h.end)) <= dist ||
      fabs(double(cluster_end) - double(h.end)) <= dist ||
      fabs(double(cluster_end) - double(h.begin)) <= dist) {
      return 1;
    }
    return 0;
  }
  return 0;
}

void PrintHit(Hit hit) {
  cout << hit.name << "\t" << hit.target << "\t" << hit.orientation <<
    "\t" << hit.begin << "\t" << hit.end << "\t" << hit.score << endl;
}
