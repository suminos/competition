///////////////////////////////////////////////////////////////////////////////
// Aho-Corasick
// キーワードが複数ある場合に，テキストの長さのみに線形な時間でパターンマッチングを行う．
// 前処理としてパターンマッチングオートマトンを作成する．(O(キーワードの数+キーワードの長さ))
// 下記コードでは，next関数やstatus関数の実装を横着しているのでワーストで数倍程度の時間がかかるかも
// 実行時間がシビアな場合は，前処理を追加することで高速化できるかもしれない．

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <queue>
#include <cstdio>
using namespace std;
#define rep(i,n) for(int i=0; i<(int)(n); i++)
#define all(c) (c).begin(), (c).end()

struct PMA {
  const static int nsize = 26;
  int c2i(char c) { return c - 'a'; }
  char i2c(int i) { return 'a' + i; }
  
  struct Node {
    int next[nsize], fail, match;
    Node() { rep(i, nsize) next[i] = 0; fail = 0; match = 0; }
  };
  vector<Node> Tree;
  
  PMA(vector<string> key) {
    Tree.push_back(Node());
    rep(i, key.size()) _add(key[i], i+1);
    _set_fail();
  }
  
  void _add(string s, int val = 1) {
    int cur = 0;
    rep(i, s.size()) {
      if(!Tree[cur].next[c2i(s[i])]) {
        Tree[cur].next[c2i(s[i])] = Tree.size();
        Tree.push_back(Node());
      }
      cur = Tree[cur].next[c2i(s[i])];
    }
    Tree[cur].match = val;
  }
  void _set_fail() {
    queue<int> qu;
    rep(i, nsize) if(Tree[0].next[i]) qu.push(Tree[0].next[i]);
    while(!qu.empty()) {
      int cur = qu.front(); qu.pop();
      rep(i, nsize) if(Tree[cur].next[i]) {
        Tree[Tree[cur].next[i]].fail = next(Tree[cur].fail, i2c(i));
        qu.push(Tree[cur].next[i]);
      }
    }
  }
  
  int next(int cur, char c) {
    if(Tree[cur].next[c2i(c)])
      return Tree[cur].next[c2i(c)];
    if(cur == 0) return 0;
    return next(Tree[cur].fail, c);
  }
  vector<int> status(int cur) {
    vector<int> r;
    while(cur) {
      if(Tree[cur].match) r.push_back(Tree[cur].match);
      cur = Tree[cur].fail;
    }
    return r;
  }
};

int main() {
  string vs[] = { "ab", "bc", "bab", "d", "abcde" };
  PMA pma(vector<string>(vs, vs+5));
  rep(i, pma.Tree.size()) {
    cout << i << ": ";
    rep(j, pma.nsize)
      if(pma.Tree[i].next[j])
        cout << pma.i2c(j) << "-" << pma.Tree[i].next[j] << " ";
    cout << "fail-" << pma.Tree[i].fail;
    if(pma.status(i).size())
      rep(j, pma.status(i).size())
        cout << " \"" << vs[pma.status(i)[j]-1] << "\"";
    cout << endl;
  }
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Rolling Hash

ll extgcd(ll a, ll b, ll &x, ll &y) {
  int g = a; x = 1; y = 0;
  if(b != 0) g = extgcd(b, a % b, y, x), y -= (a / b) * x;
  return g;
}

ll invmod(ll a, ll m) {
  ll x, y;
  if(extgcd(a, m, x, y) == 1) return (x + m) % m;
  else              return 0;
}

struct HashedSequence {
  ll val, w, x;
  static const ll m = (1LL << 60);
  HashedSequence(ll w=27): val(0), w(w), x(1) {}
  void add_R(int c) { val+=c*x; x*=w; val%=m; x%=m; }
  void del_R(int c) { x*=invmod(w,m); val-=c*x; val%=m; x%=m; }
  void add_L(int c) { val*=w; val+=c; x*=w; val%=m; x%=m; }
  void del_L(int c) { x*=invmod(w,m); val-=c; val*=invmod(w,m); val%=m; x%=m; }
};


///////////////////////////////////////////////////////////////////////////////
// 構文解析 O(N^2)

int eval(string s) {
  const int n = s.size();
  for(int i=n-1, level=0; i>=0; i--) {
    if(s[i] == '(') level++;
    if(s[i] == ')') level--;
    if(!level && s[i] == '+') return eval(s.substr(0,i)) + eval(s.substr(i+1));
    if(!level && s[i] == '-') return eval(s.substr(0,i)) - eval(s.substr(i+1));
  }
  for(i=n-1, level=0; i>=0; i--) {
    if(s[i] == '(') level++;
    if(s[i] == ')') level--;
    if(!level && s[i] == '*') return level(s.substr(0,i)) * level(s.substr(i+1));
    if(!level && s[i] == '/') return level(s.substr(0,i)) / level(s.substr(i+1));
  }
  if(s[0] == '(') return eval(s.substr(1, n-2));
  return atoi(s.c_str());
}


///////////////////////////////////////////////////////////////////////////////
// 構文解析 O(N)



///////////////////////////////////////////////////////////////////////////////
// 文字列の置換

string replace(string s, const string &from, const string &to) {
  int pos = s.find(from);
  if(pos != -1) s.replace(pos, from.size(), to);
  return s;
}

string replace_all(string s, const string &from, const string &to) {
  for(int pos = s.find(from); pos != npos; pos = s.find(from, pos)) {
    s.replace(pos, from.size(), to);
    pos += to.size();
  }
  return s;
}

