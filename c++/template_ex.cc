#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <cstdio>
#include <cmath>
using namespace std;
#define rep(i,n) for(int i=0; i<(int)(n); i++)
#define all(c) (c).begin(), (c).end()
#define pb push_back
#define mp make_pair
typedef long long ll;
typedef pair<int, int> ii;
const int inf = 1e+9;

template <class T> string t2s(const T& t) { ostringstream os; os << t; return os.str(); }
template <class T> T s2t(const string& s) { istringstream is(s); T t; is >> t; return t; }

template <typename T> ostream& operator<<(ostream& os, const vector<T>& cs) {
  os << "[";
  rep(i, cs.size()) os << cs[i] << ", ";
  os << "]";
  return os;
}

template <typename T> ostream& operator<<(ostream& os, const vector<vector<T> >& cs) {
  int len = 0;
  rep(i, cs.size()) rep(j, cs[i].size()) len = max(len, (int)t2s(cs[i][j]).size());
  os << endl;
  rep(i, cs.size()) {
    os << "  [";
    rep(j, cs[i].size()) { os.width(len); os << cs[i][j] << ", "; }
    os << "]" << endl;
  }
  return os;
}

#define trace(x) (cerr << "\033[31m" << #x << ": " << (x) << "\033[39m" << endl)


int main() {
}
