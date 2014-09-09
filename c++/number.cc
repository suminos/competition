// エラトステネスの篩
vector<int> p;
void sieve(int n) {
  p = vector<int>(n);
  rep(i, n) p[i] = i;
  p[0] = p[1] = 0;
  for(int i=2; i*i<=n; i++)
    if(p[i])
      for(int j=i*2; j<n; j+=i)
        p[j] = 0;
  p.erase(remove(all(p), 0), p.end());
}

// ユークリッドの互除法 O(log n)
int gcd(int a, int b) {
  return b != 0 ? gcd(b, a%b) : a;
}

int lcm(int a, int b) {
  return a * b / gcd(a, b);
}

/*
  gcd は #include <algorithm> して， __gcd(a,b) を使う．
  lcm は a/__gcd(a,b)*b とする．
*/


// ax + by = gcd(a, b)
int extgcd(int a, int b, int& x, int& y) {
  int g = a; x = 1; y = 0;
  if(b != 0) g = extgcd(b, a%b, y, x), y -= (a / b) * x;
  return g;
}

// オイラーのφ関数 (verified: PKU2407)
// nと互いに素な整数の個数を返す．
int euler(int n) {
  if(n == 0) return 0;
  int r = n;
  for(int x=2; x*x<=n; x++) {
    if(n % x == 0) {
      r -= r / x;
      while(n % x == 0) n /= x;
    }
  }
  if(n > 1) r -= r / n;
  return r;
}

// ax + z = y となる a,z
// ここで x,y \in Q, a \in N, z \in Q, 0<=z<abs(x)
double mod(double x, double y) {
  int a = int(y / x);
  double z = y - a * x;
  return z;
}

// 総和を取るとオーバーフローする場合の平均のとり方
ll average(vector<ll> xs) {
  ll t = 0, avg = 0;
  rep(i, n) {
    t += xs[i];
    avg += t / n;
    t %= n;
  }
  return avg;
}

// 組み合わせをDPで求める
vector<vector<int> > C;
void init(int n) {
  C = vector<vector<int> >(n, vector<int>(n));
  rep(i,n) C[i][0] = C[i][i] = 1;
  rep(i,n) for(int j=1; j<i; j++) C[i][j] = (C[i-1][j-1]+C[i-1][j]) % mod;
}

// 三分探索
double tridirectional_search(double (*f) (double)) {
  double lb = -inf, ub = inf;
  rep(i, 100) {
    double m1 = lb + (ub-lb) / 3.0;
    double m2 = lb + (ub-lb) / 3.0 * 2.0;
    if(f(m1) > f(m2)) ub = m2;
    else lb = m1;
  }
  return lb < 0.0 ? 0.0 : lb > 1.0 ? 1.0 : lb;
}

// Xorshift - 擬似乱数生成アルゴリズム(http://ja.wikipedia.org/wiki/Xorshift)
typedef unsigned int uint;
uint xor128() {
  static uint x = 123456789;
  static uint y = 362436069;
  static uint z = 521288629;
  static uint w = 88675123; 
  uint t;
  
  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)); 
}

// 組み合わせ列挙 - (http://d.hatena.ne.jp/jetbead/20121202/1354406422)
// (n,k)=(4,2)のとき，0011,0101,0110,1001,1010,1100
void print_bit(ull S, int n=64){
  for(int i=n-1; i>=0; i--){
    if(S>>i & 1) std::cout << 1;
    else std::cout << 0;
  }
  std::cout << std::endl;
}
 
void subset_combination(int n, int k){
  ull S = (1ULL << k) - 1ULL;
  ull E = ~((1ULL << n) - 1ULL);
  while(!(S & E)){
    print_bit(S, n);
    ull smallest = S & -S;
    ull ripple = S + smallest;
    ull nsmallest = ripple & -ripple;
    S = ripple | (((nsmallest / smallest) >> 1) - 1);
  }
}
 
int main(){
  subset_combination(4, 2);
  return 0;
}

// ガウスの除去法
// ( ax + by ) = (c)
// ( dx + ey ) = (f) の形式で，Nx(N+1)行列を与える
vector<double> gaussian_elimination(vector<vector<double> > A) {
  const int n = A.size();
  vector<double> r(n);
  
  rep(i, n) {
    {
      double t = A[i][i];
      rep(j, n+1) A[i][j] /= t;
    }
    for(int k=i+1; k<n; k++) {
      double t = A[k][i];
      rep(j, n+1) A[k][j] -= t * A[i][j];
    }
  }
  
  for(int i=n-1; i>=0; i--) {
    for(int j=i+1; j<n; j++) A[i][n] -= A[i][j] * r[j];
    r[i] = A[i][n];
  }
  
  return r;
}

