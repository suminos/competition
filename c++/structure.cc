struct UnionFind {
  vector<int> v;
  UnionFind(int n): v(n) {
    rep(i, n) v[i] = i;
  }
  int find(int x) {
    if(v[x] == x) return x;
    return v[x] = find(v[x]);
  }
  void unite(int x, int y) {
    v[find(x)] = find(y);
  }
  bool same(int x, int y) {
    return find(x) == find(y);
  }
};

struct SegTree {
  vector<int> v;
  SegTree(int n): v(n*2-1, inf) {} // nは2の累乗であること
  void update(int i, int x) {
    i += n - 1;
    v[i] = x;
    while(i > 0) {
      i = (i - 1) / 2;
      v[i] = min(v[i * 2 + 1], v[i * 2 + 2]);
    }
  }
  // 外からは query(a, b, 0, 0, n) として呼ぶ．
  int query(int a, int b, int k, int l, int r) {
    if(r <= a || b <= l) return inf;
    if(a <= l && r <= b) return v[k];
    return min( query(a, b, k * 2 + 1, l, (l + r) / 2),
                query(a, b, k * 2 + 2, (l + r) / 2, r) );
  }
};


// i \in [1,n]
struct FenwickTree {
  int n;
  vector<int> v;
  
  FenwickTree(int n): n(n), v(n+1) {}
  int sum(int i) {
    int s = 0;
    while(i > 0) {
      s += v[i];
      i -= i & -i;
    }
    return s;
  }
  
  void add(int i, int x) {
    while(i <= n) {
      v[i] += x;
      i += i & -i;
    }
  }
};

// unverified
// x \in [1,w], y \in [1,h]
struct FenwickTree2D {
  int n;
  vector<FenwickTree> v;
  
  FenwickTree(int w, int h): w(w), h(h), v(h+1, vector<FenwickTree>(w)) {}
  int sum(int x, int y) {
    int s = 0;
    while(y > 0) {
      s += v[y].sum(x);
      y -= y & -y;
    }
    return s;
  }
  
  void add(int x, int y, int val) {
    while(y <= n) {
      v[y].add(x) += val;
      y += y & -y;
    }
};

struct RollingHash {
};

struct D_Tree {
};

