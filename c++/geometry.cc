#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;
#define rep(i,n) for(int i=0; i<(int)(n); i++)
#define all(c) (c).begin(), (c).end()
#define iter(c) __typeof((c).begin())
typedef long long ll;
const int inf = (1<<28);

const double eps = 1e-8;
const double pi = acos(-1.0);

typedef complex<double> P;

struct L {
  P fst, snd;
  L() {}
  L(P fst, P snd) : fst(fst), snd(snd) {}
};

struct LS {
  P fst, snd;
  LS() {}
  LS(P fst, P snd) : fst(fst), snd(snd) {}
};

struct C {
  P p;
  double r;
  C() {}
  C(P p, double r) : p(p), r(r) {}
};

// 内積 |a||b|cos
double dot(P a, P b) { return real( conj(a) * b ); }

// 外積 |a||b|sin
double cross(P a, P b) { return imag( conj(a) * b ); }

// 回転方向関数 a->b からの c の方向
int ccw(P a, P b, P c) {
  b -= a; c -= a;
  double v = cross(b, c);
  if (0.0 < v)           return  1; // c-clockwise
  if (v < 0.0)           return -1; // clockwise
  if (dot(b,c) < 0.0)    return  2; // c--a--b
  if (norm(b) < norm(c)) return -2; // a--b--c
  return 0;                        // a--c--b
}

P middle_point(P a, P b, double m = 0.5) {
  return a + (b - a) * m;
}

double dist(P a, L bc) {
  P b = bc.fst, c = bc.snd;
  return abs( cross(a-b, c-b) ) / abs(c - b);
}

double dist(P a, LS bc) {
  P b = bc.fst, c = bc.snd;
  if (dot(a-b, c-b) < eps) return abs(a - b);
  if (dot(a-c, b-c) < eps) return abs(a - c);
  return abs( cross(a-b, c-b) ) / abs(c - b);
}

double dist(L ab, L cd) {
  P a = ab.fst, b = ab.snd, c = cd.fst, d = cd.snd;
  if (abs( cross(b-a, d-c) ) > 0.0) return 0.0;
  return abs( cross(a-c, d-c) ) / abs(d - c);
}

// require : dist(P,L)
double dist(L ab, LS cd) {
  P a = ab.fst, b = ab.snd, c = cd.fst, d = cd.snd;
  if (cross(b-a, c-a) * cross(b-a, d-a) < 0.0) return 0.0;
  return min( dist(c, ab), dist(d, ab) );
}

double dist(L ab, C c) {
  double d = dist(ab, c.p);
  return d < c.r ? 0.0 : d;
}

// require : dist(P,LS)
double dist(LS ab, LS cd) {
  P a = ab.fst, b = ab.snd, c = cd.fst, d = cd.snd;
  if ( cross(b-a, c-a) * cross(b-a, d-a) < 0.0 &&
       cross(d-c, a-c) * cross(d-c, b-c) < 0.0 ) return 0.0;
  double res = dist(a, cd);
  res = min( res, dist(b, cd) );
  res = min( res, dist(c, ab) );
  res = min( res, dist(d, ab) );
  return res;
}

double dist(LS ab, C c) {
  double abc = dist(c.p, ab);
  double ac = abs(c.p - ab.fst);
  double bc = abs(c.p - ab.snd);
  if (ac > bc) swap(ac,bc);
  if (c.r < abc)            return abc - c.r; // 線分が円の外にある
  if (ac < c.r && c.r < bc) return 0.0;       // 一点で交わる
  if (bc < c.r)             return c.r - bc;  // 線分が円の中にある
  return 0.0;                                 // 二点で交わる
}

P cross_point(L ab, L cd) {
  P a = ab.fst, b = ab.snd, c = cd.fst, d = cd.snd;
  return a + ( cross(d-c, c-a) / cross(d-c, b-a) ) * (b - a);
}

P cross_point(L ab, LS cd) {
  P a = ab.fst, b = ab.snd, c = cd.fst, d = cd.snd;
  return a + ( cross(d-c, c-a) / cross(d-c, b-a) ) * (b - a);
}

P cross_point(LS ab, LS cd) {
  P a = ab.fst, b = ab.snd, c = cd.fst, d = cd.snd;
  return a + ( cross(d-c, c-a) / cross(d-c, b-a) ) * (b - a);
}

pair<P,P> cross_point(L ab, C c) {
  P a = ab.fst, b = ab.snd;
  P p = P(b-a) / abs(b-a) * polar(1.0, pi/2.0);
  double t = acos(dist(c.p, ab) / c.r);
  return pair<P,P>( c.p + p * polar(c.r,  t),
                    c.p + p * polar(c.t, -t) );
}

// require : cross_point(L,C), dist(P,LS)
pair<P,P> cross_point(LS ab, C c) {
  pair<P,P> pp = intersect( L(ab.fst, ab.snd), c );
  if ( dist(pp.first , ab) == 0.0 ) pp.first  = P(nan,nan);
  if ( dist(pp.second, ab) == 0.0 ) pp.second = P(nan,nan);
  return pp;
}

pair<P,P> cross_point(C a, C b) {
  P c1 = a.p, c2 = b.p;
  double r1 = a.r, r2 = b.r;
  double d = abs(c1 - c2);
  double s = (d*d + r1*r1 - r2*r2) / (2.0*d);
  double t = sqrt(r1*r1 - s*s);
  P diff = (c2 - c1) / d;
  return pair<P,P> ( c1 + diff * P(s,  t),
                     c1 + diff * P(s, -t) );
}

// 符号付き面積がほしい時は abs() を取る
double area(vector<P> ps) {
  double a = 0.0;
  rep (i, ps.size()) a += cross( ps[i], ps[(i+1)%ps.size()] );
  return abs(a) / 2.0;
}

// 円と円の共通部分の面積
double area(const C& c1, const C& c2) {
    double d = abs(c1.p - c2.p);
    if (c1.r + c2.r <= d + eps) {
        return 0.0;
    } else if (d <= abs(c1.r - c2.r) + eps) {
        double r = min(c1.r, c2.r);
        return r * r * pi;
    } else {
        double rc = (d*d + c1.r*c1.r - c2.r*c2.r) / (2*d);
        double tha = acos( rc      / c1.r);
        double thb = acos((d - rc) / c2.r);
        return c1.r*c1.r*tha + c2.r*c2.r*thb - d*c1.r*sin(tha);
    }
}

// require : dist(P,LS)
bool is_on(P a, vector<P> ps) {
  rep (i, ps.size())
    if ( dist(a, LS( ps[i], ps[(i+1)%ps.size()]) ) < eps ) return true;
  return false;
}

// ps は凸多角形であること．（実際にはaから見て凸多角形なら正しく判定できる）
bool is_in(P a, vector<P> ps) {
  rep (i, ps.size()) ps[i] -= a;
  double sign = cross(ps[0], ps[1]);
  rep (i, ps.size())
    if (sign * cross( ps[i], ps[(i+1)%ps.size()] ) < 0.0)
      return false;
  return true;
}

// require : dist(LS,LS)
// ps が非凸でも，ねじれていてもOK
// is_on 判定はしておくべき
bool is_in(P a, vector<P> ps) {
  double r = ((double)rand() / RAND_MAX) * 2.0 * pi; // 乱数でごまかす
  P b = a + polar((double)inf, r);
  int cnt = 0;
  rep (i, ps.size()) {
    P c = ps[i], d = ps[(i+1) % ps.size()];
    if (dist( LS(a,b), LS(c,d) ) < eps) cnt++;
  }
  return cnt % 2 == 1;
}

// usage : sort(all(ps), xycomp)
// (x座標, y座標) でソート
bool xycomp(const P &a, const P &b) {
  if (a.real() != b.real()) return a.real() < b.real();
  return a.imag() < b.imag();
}
// (偏角, 長さ) でソート
bool argcomp(const P &a, const P&b) {
  double aa = arg(a) < 0.0 ? 2*pi-arg(a) : arg(a);
  double bb = arg(b) < 0.0 ? 2*pi-arg(b) : arg(b);
  if(aa != bb) return aa < bb;
  return abs(a) < abs(b);
}

// Andrew scan : O(n log n)
// require : ccw, xycomp
// 凸包の辺上の点も含めたい場合は， ccw()<0 を ccw()==-1 に書き換える
vector<P> convex_hull(vector<P> ps) {
  int n = ps.size(), k = 0;
  sort(all(ps), xycomp);
  vector<P> ch(2*n);
  for (int i = 0; i < n; ch[k++] = ps[i++])              // lower
    while (k>=2 && ccw(ch[k-2], ch[k-1], ps[i]) < 0) --k;
  for (int i = n-2, t = k+1; i >= 0; ch[k++] = ps[i--]) // upper
    while (k >= t && ccw(ch[k-2], ch[k-1], ps[i]) < 0) --k;
  ch.resize(k-1);
  return ch;
}

pair<L,L> tangent(P a, C c) {
  double t = asin(c.r / abs(c.p - a));
  return pair<L,L> ( L( a, (c.p - a) * polar(1.0,  t) ),
                     L( a, (c.p - a) * polar(1.0, -t) ) );
}

// ( 内接線のペア, 外接線のペア ) を返す
// a1 や b1 は接点．二円は交差したり接したりしないこと．
// ランダム入力でテスト済み
pair< pair<L,L>, pair<L,L> > tangents(C a, C b) {
  const P c1 = a.p, c2 = b.p;
  const double r1 = a.r, r2 = b.r;
  const P p = (c2 - c1) / abs(c2 - c1);
  
  double t1 = acos((r2 - r1) / abs(c2 - c1));
  P a1 = c1 + p * polar(r1, pi-t1);
  P b1 = c2 + p * polar(r2, pi-t1);
  P a2 = c1 + p * polar(r1, pi+t1);
  P b2 = c2 + p * polar(r2, pi+t1);
  pair<L,L> ll1 = pair<L,L> ( L(a1, b1), L(a2, b2) );
  
  double t2 = acos((r1 + r2) / abs(c2 - c1));
  P a3 = c1 + p * polar(r1,  t2);
  P b3 = c2 - p * polar(r2,  t2);
  P a4 = c1 + p * polar(r1, -t2);
  P b4 = c2 - p * polar(r2, -t2);
  pair<L,L> ll2 = pair<L,L> ( L(a3, b3), L(a4, b4) );
  
  return pair< pair<L,L>, pair<L,L> > (ll1, ll2);
}

// 三角形の外心
// 辺の垂直二等分線を二つ求め，交点を返す
// require: P, L, crosspoint(L, L)
P circumcenter(P a, P b, P c) {
  P x = (a + b) / 2.0, y = (a + c) / 2.0;
  L l1 = L( x, x + (b-a) * polar(1.0, pi / 2.0) );
  L l2 = L( y, y + (c-a) * polar(1.0, pi / 2.0) );
  return cross_point(l1, l2);
}

// 三角形の外接円
// 辺の垂直二等分線を二つ求め，交点を返す
// require: P, L, C, crosspoint(L, L)
C circumcircle(P a, P b, P c) {
  P x = (a + b) / 2.0, y = (a + c) / 2.0;
  L l1 = L( x, x + (b-a) * polar(1.0, pi / 2.0) );
  L l2 = L( y, y + (c-a) * polar(1.0, pi / 2.0) );
  P cent = cross_point(l1, l2);
  return C( cent, abs(cent - a) );
}

// 多角形の外接円
// 三角形の外接円を求める場合もこちらでOK
C circumcircle(const vector<P> &ps) {
  C c = C( P(0.0, 0.0), 1.0 );
  double move = 0.5;
  while( move > eps ) {
    rep(k, 100) {
      int x = 0;
      rep(i, ps.size())
        if(norm(ps[i] - c.p) > norm(ps[x] - c.p)) x = i;
      c.p += (ps[x] - c.p) * move;
      c.r = abs(ps[x] - c.p);
    }
    move /= 2.0;
  }
  return c;
}

/*
a,b,cは基本的に辺の長さ．

  ヘロンの公式
  S = sqrt(s(s-a)(s-b)(s-c))  (s=三角形の周の長さの半分)
  
  s = abc/4R = 2r/(a+b+c)     (R=外接円の半径，r=内接円の半径)
  
  a^2 = b^2+c^2-2bc*cos(A)    (a,b,cについて対称)
  cos(A) = (b^2+c^2-a^2)/2bc  (上の式の変形)
  
  R = a/2sin(A)               
    = abc/2cross(b-a,c-a)     (crossに渡しているa,b,cはベクトル)

  ピックの定理
    S = i + b/2 - 1           (S=多角形の面積， i=内部格子点, b=辺上格子点)
 */

////////////////////////////////////////////
// 空間幾何
// valarray を使う．
// http://www.prefield.com/algorithm/geometry3D/geometry3D.html

typedef valarray<double> P;

double dot(const P &a, const P &b) {
  return (a * b).sum();
}

P cross(const P &a, const P &b) {
  return a.cshift(+1) * b.cshift(-1) -
          a.cshift(-1) * b.cshift(+1);
}

double dist(const P &a, const P &b) {
  return dot(a-b, a-b);
}

