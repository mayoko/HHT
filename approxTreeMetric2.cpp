// 距離空間が与えられるので, 近似的に木距離にするプログラム
// O(N^2)

#include<iostream>
#include<random>
#include<vector>
#include<algorithm>
#include<unordered_map>
#include<cassert>
#include<cstring>

using namespace std;

random_device rnd;
mt19937 mt(rnd());
// 距離の最大値
const int MAX = 114514;
// 頂点数の最大値
const int MAXV = 5555;
// 頂点数(入力)
int N;

struct Edge {
    int to;
    int cost;
    Edge() {}
    Edge(int t, double c) : to(t), cost(c) {}
};

// 距離空間
int dist[MAXV][MAXV];
// 距離空間で最大の距離と最小の距離
int mini = MAX, maxi = 0;
// 木が何段組みになるかってヤツ
int delta;
// 例の乱数
double beta;
// 例のランダム順列
int perm[MAXV];

// 木(出力すべきもの)
vector<Edge> T[2*MAXV];
// 木の頂点に対応する頂点群がシングルトンになっているとき, 対応する頂点をメモしておく(これも出力すべきもの)
int corr[2*MAXV];
// ある頂点が, 今順列の何番目まで進んでいるか
int cur[MAXV];
// 今のところの木の頂点数
int treeNum;

// 深さ depth のところで頂点群 vs が残っているので適当に分割して木を作る
// 集合 vs に対応する木の頂点は node とラベリングされてる
void dfs(const vector<int>& vs, int depth, int node) {
    // 新しい頂点についたので
    ++treeNum;
    if (vs.size() == 1) {
        corr[node] = vs[0];
        return;
    }
    // vs を最も近い頂点でグループ分け
    unordered_map<int, vector<int> > mp;
    double radius = beta * pow(2, delta-depth-1)*mini;
    int cost = pow(2, delta-depth)*mini;
    for (int el : vs) {
        int now = cur[el];
        while (now < N) {
            if (dist[el][perm[now]] <= radius) break;
            ++now;
        }
        assert(now < N);
        mp[perm[now]].push_back(el);
    }
    for (auto p : mp) {
        T[node].emplace_back(treeNum, cost);
        dfs(p.second, depth+1, treeNum);
    }
}

int main() {
    cin >> N;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> dist[i][j];
        }
    }
    // 最小値, 最大値を持ってくる
    mini = MAX, maxi = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) mini = min(mini, dist[i][j]);
            maxi = max(maxi, dist[i][j]);
        }
    }
    // 木が何段必要なのかを求める
    {
        int now = mini;
        while (1) {
            if (now >= maxi) break;
            now *= 2;
            delta++;
        }
    }
    cerr << delta << endl;
    // ランダム順列と乱数を用意しておく
    for (int i = 0; i < N; i++)
        perm[i] = i;
    shuffle(perm, perm+N, mt19937());
    // とりあえず一様乱数で良いや
    // TODO: 1/(x ln 2) の形の乱数を用意する
    uniform_real_distribution<double> uni(1.0, 2.0);
    beta = uni(mt);
    // 木の作成
    vector<int> vs(N);
    for (int i = 0; i < N; i++)
        vs[i] = i;
    memset(corr, -1, sizeof(corr));
    dfs(vs, 0, 0);
    // 木の出力
    cout << treeNum << endl;
    for (int i = 0; i < treeNum; i++) {
        cout << i << " " << corr[i] << endl;
        cout << T[i].size() << endl;
        for (Edge e : T[i]) {
            cout << e.to << " " << e.cost << endl;
        }
        cout << endl;
    }
}