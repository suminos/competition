///////////////////////////////////////////////////////////////////////////////
// トポロジカルソート

void _dfs(const Graph &g, int cur, vecto<int> &ord, int &num) {
  rep(i, g[cur].size()) {
    int nxt = g[cur][i];
    if(ord[nxt] != -1) continue;
    _dfs(g, nxt, ord, num);
  }
  ord[cur] = --num;
}

// [ トポロジカルソート ] : O(V+E)
// input  : DAG
// output : 辺 (u,v) が存在するなら o[u]<o[v] となる長さ V の配列
vector<int> topological_sort(const Graph &g) {
  const int n = g.size();
  
  int num = n;
  vector<int> ord(n, -1);
  
  rep(i, n) if(ord[i] == -1) _dfs(g, i, ord, num);
  
  vector<int> _ord(n);
  rep(i, n) _ord[ord[i]] = i;
  return _ord;
}

///////////////////////////////////////////////////////////////////////////////
// 二部グラフの最大マッチング
// required: 
// verified: X

bool _dfs(const Graph &g, int n, int m, int u, vector<int> &X, vector<int> &Y, vector<int> &used) {
  if(used[u]) return false;
  used[u] = true;
  rep(k, g[u].size()) {
    int v = g[u][k] - n;
    if(Y[v] == -1) { X[u] = v; Y[v] = u; return true; }
    else if(_dfs(g, n, m, Y[v], X, Y, used)) { X[u] = v; Y[v] = u; return true; }
  }
  return false;
}

// 二部グラフ K_n,m の最大マッチングを求める．
// g.size() == n + m であること．
int bipartite_matching(const Graph &g, int n, int m, vector<int> &X, vector<int> &Y) {
  X = vector<int>(n, -1), Y = vector<int>(m, -1); // 分割の要素のマッチング相手
  rep(i, n) {
    vector<int> used(n, false);
    _dfs(g, n, m, i, X, Y, used);
  }
  return n - count(all(X), -1);
}

///////////////////////////////////////////////////////////////////////////////
// 最大流
// required: <queue>, inf

struct Edge {
  int dst, res;
  Edge() {}
  Edge(int dst, int res): dst(dst), res(res) {}
};

typedef vector<vector<int> > Graph;
typedef vector<Edge> Edges;

void add_edge(Graph &g, Edges &es, int src, int dst, int capgo, int capbk=0) {
  g[src].push_back(es.size());
  es.push_back(Edge(dst, capgo));
  g[dst].push_back(es.size());
  es.push_back(Edge(src, capbk));
}

// levelize residue graph
bool bfs(const Graph &g, const Edges &es, int s, int t, vector<int> &level) {
  rep(i, g.size()) level[i] = -1;
  level[s] = 0;
  queue<int> qu;
  qu.push(s);
  while(!qu.empty()) {
    int cur = qu.front(); qu.pop();
    rep(i, g[cur].size()) {
      if(level[es[g[cur][i]].dst] != -1) continue;
      if(es[g[cur][i]].res <= 0) continue;
      level[es[g[cur][i]].dst] = level[cur] + 1;
      qu.push(es[g[cur][i]].dst);
    }
  }
  return level[t]!=-1;
}

// find augument path
bool dfs(const Graph &g, const Edges &es, const vector<int> &level, int s, int t, vector<int> &path) {
  const int n = g.size();
  if(s == t) return true;
  rep(i, g[s].size()) {
    if(level[es[g[s][i]].dst] != level[s]+1) continue;
    if(es[g[s][i]].res == 0) continue;
    path.push_back(g[s][i]);
    if(dfs(g, es, level, es[g[s][i]].dst, t, path)) return true;
    path.pop_back();
  }
  return false;
}

// Dinic 法の隣接リスト実装 : O(V^2 * E)
// es[i] の逆辺を es[i^1] で表す．
int max_flow(const Graph &g, Edges es, int s, int t) {
  const int n = g.size();
  int total = 0;
  vector<int> level(n), path;
  while(bfs(g, es, s, t, level) && dfs(g, es, level, s, t, path)) {
    int aug = inf;
    rep(i, path.size()) aug = min(aug, es[path[i]].res);
    rep(i, path.size()) es[path[i]].res -= aug;
    rep(i, path.size()) es[path[i]^1].res += aug;
    total += aug;
    path.clear();
  }
  return total;
}

///////////////////////////////////////////////////////////////////////////////
// 頂点連結度
// required: max_flow

int vertex_connectivity(const Graph &g, int src, int dst) {
  const int n = g.size();
  vector<int> ass(n);
  for(int i=0,j=0; i<n; i++) {
    if(i==src || i==dst) ass[i] = -1;
    else ass[i] = j++;
  }
  Graph gg((n-2)*2+2);
  Edges es;
  rep(i,n) if(!(i==src || i==dst))
    add_edge(gg,es,ass[i]*2+2,ass[i]*2+2+1,1);
  rep(i,n) rep(jj,g[i].size()) {
    int j = g[i][jj];
    int s = i==src ? 0 : i==dst ? 1 : ass[i]*2+2+1;
    int t = j==src ? 0 : j==dst ? 1 : ass[j]*2+2;
    add_edge(gg,es,s,t,1);
  }
  
  return max_flow(gg,es,0,1);
}

int vertex_connectivity(const Graph &g) {
  const int n = g.size();
  int v = 0;
  rep(i,n) if(g[i].size()<g[v].size()) v = i;
  vector<vector<bool> > adj(n,vector<bool>(n,false));
  rep(i,n) rep(jj,g[i].size()) adj[i][g[i][jj]] = true;
  int k = inf;
  rep(i,n) if(i!=v && !adj[v][i])
    k = min(k, vertex_connectivity(g,v,i));
  rep(i,g[v].size()) rep(j,i) if(!adj[g[v][i]][g[v][j]])
    k = min(k, vertex_connectivity(g,g[v][i],g[v][j]));
  return k;
}

///////////////////////////////////////////////////////////////////////////////
// 有向オイラー路 : O(n+m)
// 先に連結判定の必要なし
// path には頂点列が入る

typedef vector<vector<int> > Graph;

void visit(Graph& g, int u, vector<int>& path) {
  while(!g[u].empty()){
    int v = g[u].back();
    g[u].pop_back();
    visit(g, v, path);
  }
  path.push_back(u);
}

bool euler_path(g, int s, vector<int> &path) {
  int n = g.size(), m = 0;
  vector<int> deg(n, 0); // 出次数を数える
  rep(i, n){
    m += g[i].size();
    deg[i] += g[i].size();       // out
    rep(j, g[i]) deg[g[i][j]]--; // in
  }
  int k = n - count(all(deg), 0);
  if (k == 0 || (k == 2 && deg[s] == 1)) {
    path.clear();
    dfs(g, s, path);
    reverse(all(path));
    return path.size() == m + 1;
  }
  return false;
}

///////////////////////////////////////////////////////////////////////////////
// 無向オイラー路 : O(n^2)
// 内部で隣接リストから隣接行列への変換を行なっているため

typedef vector<vector<int> > Graph;

void visit(const Graph &g, vector<vector<int> > &adj, int u, vector<int> &path) {
  rep(i, g[u]) if(adj[u][g[u][i]]) {
    adj[u][g[u][i]]--;
    adj[g[u][i]][u]--;
    visit(g, adj, g[u][i], path);
  }
  path.push_back(u);
}

bool euler_path(Graph g, int s, vector<int> &path) {
  int n = g.size(), m = 0, odd = 0;
  REP(i, n) {
    if(g[i].size() % 2 == 1) odd++;
    m += g[i].size();
  }
  m /= 2;
  if(odd == 0 || (odd == 2 && g[s].size() % 2 == 1)) {
    vector<vector<int> > adj(n, vector<int>(n, 0));
    rep(i, n) rep(j, g[i]) adj[i][g[i][j]];
    path.clear();
    visit(g, adj, s, path);
    reverse(all(path));
    return path.size() == m + 1;
  }
  return false;
}

///////////////////////////////////////////////////////////////////////////////
// 強連結成分分解 計算量:O(V+E) verified:PKU1236?

void _dfs(const Graph &g, int cur, vector<int> &num, vector<int> &low, stack<int> &st, int &o, vector<int> &scc) {
  num[cur] = low[cur] = o++;
  st.push(cur);
  rep(i, g[cur].size()) {
    int nxt = g[cur][i];
    if(num[nxt] == -1) {
      _dfs(g, nxt, num, low, st, o, scc);
      low[cur] = min(low[cur], low[nxt]);
    }
    else if(scc[nxt] == -1) {
      low[cur] = min(low[cur], low[nxt]);
    }
  }
  if(low[cur] == num[cur]) {
    while(true) {
      int t = st.top(); st.pop();
      scc[t] = cur;
      if(t == cur) break;
    }
  }
}

vector<int> scc_decomp(const Graph &g) {
  const int n = g.size();
  vector<int> scc(n, -1);
  vector<int> num(n, -1), low(n);
  stack<int> st;
  int o = 0;
  rep(i, n) if(num[i] == -1) _dfs(g, i, num, low, st, o, scc);
  return scc;
}

///////////////////////////////////////////////////////////////////////////////
// 最小費用流 計算量:O(V^2 F) verified:AOJ:2429
// グラフは有向グラフであること．無向グラフの場合は頂点を倍にして対応する
// 任意の辺(src, dst, cap, cst)に対して逆辺(dst,src,0,-cst)を張っておくこと．
// 返却値は(実際に流せた流量,かかったコスト)

// バグがあるので使わないこと

#define add_edge(s,t,cap,cst) (capa[s][t]+=cap, cost[s][t]+=cst, cost[t][s]-=cst)

pair<int, int> min_cost_flow(Mat cap, Mat cost, int s, int t, int f) {
  int n = cap.size();
  vector<int> pot(n, 0);
  int tf = 0, tc = 0;
  while(tf < f) {
    vector<int> dist(n, inf);
    vector<int> prev(n, -1);
    vector<bool> visited(n, false);
    dist[s] = 0;
    while(true) {
      int cur = -1;
      rep(i, n) if(!visited[i] && (cur == -1 || dist[i] < dist[cur])) cur = i;
      if(cur == -1) break;
      visited[cur] = true;
      rep(nxt, n) if(!visited[nxt] && cap[cur][nxt] > 0) {
        if(dist[nxt] > dist[cur] + cost[cur][nxt] + pot[cur] - pot[nxt]) {
          dist[nxt] = dist[cur] + cost[cur][nxt] + pot[cur] - pot[nxt];
          prev[nxt] = cur;
        }
      }
    }
    if(!visited[t]) break;
    int ff = f - tf;
    for(int p=t; p!=s; p=prev[p]) ff = min(ff, cap[prev[p]][p]);
    for(int p=t; p!=s; p=prev[p]) cap[prev[p]][p] -= ff;
    for(int p=t; p!=s; p=prev[p]) cap[p][prev[p]] += ff;
    for(int p=t; p!=s; p=prev[p]) tc += ff * cost[prev[p]][p];
    tf += ff;
    rep(i, n) pot[i] += dist[i];
  }
  return pair<int, int>(tf, tc);
}

///////////////////////////////////////////////////////////////////////////////
// 最小費用流 : O(nmf)
// 最短経路の計算にベルマンフォードを用いている
// 負閉路の検出もする場合は追記が必要

struct Edge {
  int src, dst, cap, cst;
  Edge() {}
  Edge(int src, int dst, int cap, int cst): src(src), dst(dst), cap(cap), cst(cst) {}
};

const int n_max = 1000;
typedef vector<Edge> Graph;

void add_edge(Graph &g, int src, int dst, int cap, int cst) {
  int n = g.size();
  g.pb(Edge(src, dst, cap,  cst));
  g.pb(Edge(dst, src,   0, -cst));
}

// 最小費用流: O(nmf)
ii min_cost_flow(Graph g, int src, int dst, int flow) {
  int f = 0, c = 0;
  
  while(f < flow) {
    vector<int> dist(n_max, inf); dist[src] = 0;
    vector<int> hist(n_max);

    while(true) {
      bool updated = false;
      rep(i, g.size()) {
        if(g[i].cap > 0 && dist[g[i].dst] > dist[g[i].src] + g[i].cst) {
          dist[g[i].dst] = dist[g[i].src] + g[i].cst;
          hist[g[i].dst] = i;
          updated = true;
        }
      }
      if(!updated) break;
    }
    
    if(dist[dst] == inf) return ii(f, c);
    
    int ff = inf;
    for(int p=dst; p!=src; p=g[hist[p]].src) ff = min(ff, g[hist[p]].cap);
    for(int p=dst; p!=src; p=g[hist[p]].src) g[hist[p]  ].cap -= ff;
    for(int p=dst; p!=src; p=g[hist[p]].src) g[hist[p]^1].cap += ff;
    for(int p=dst; p!=src; p=g[hist[p]].src) c += ff * g[hist[p]].cst;
    f += ff;
  }

  return ii(f, c);
}

///////////////////////////////////////////////////////////////////////////////
// 閉路存在判定
bool dfs(int cur, const Graph &g, vector<int> &col, int c) {
  col[cur] = c;
  rep(i, g[cur].size()) {
    int nxt = g[cur][i];
    if(col[nxt] == c) return true;
    if(dfs(nxt, g, col, c)) return true;
  }
  col[cur] = -c;
  return false;
}

bool is_cyclic(const Graph &g) {
  const int n = g.size();
  vector<int> col(n, inf);
  rep(i, n) if(col[i] == inf) if(dfs(i, g, col, i+1)) return true;
  return false;
}


