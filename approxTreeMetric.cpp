// 距離空間が与えられるので, 近似的に木距離にするプログラム
// O(N^2)(らしい)
// TODO: 今書いてるのは O(N^3) なので O(N^2) にしたい

#include<iostream>
#include<random>
#include<vector>
#include<algorithm>

using namespace std;

random_device rnd;
mt19937 mt(rnd());
const int MAX = 114514;

struct Edge {
    int to;
    double cost;
    Edge() {}
    Edge(int t, double c) : to(t), cost(c) {}
};

int main() {
    int N;
    cin >> N;
    vector<vector<int> > d(N, vector<int>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> d[i][j];
        }
    }
    // 最小値, 最大値を持ってくる
    int mini = MAX, maxi = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) mini = min(mini, d[i][j]);
            maxi = max(maxi, d[i][j]);
        }
    }
    // 木が何段必要なのかを求める
    int delta = 0;
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
    vector<int> perm(N);
    for (int i = 0; i < N; i++)
        perm[i] = i;
    shuffle(perm.begin(), perm.end(), mt19937());
    // とりあえず一様乱数で良いや
    // TODO: 1/(x ln 2) の形の乱数を用意する
    uniform_real_distribution<double> uni(1.0, 2.0);
    double beta = uni(mt);
    // 初期化
    // 頂点集合 V 全体を D[0] に入れる
    vector<vector<int> > D(1);
    for (int i = 0; i < N; i++)
        D[0].push_back(i);
    // 木の頂点数
    int node = 1;
    // 木グラフ
    vector<vector<Edge> > T(2*N);
    // 木グラフの対応する頂点(singleton になったところだけ)
    vector<int> cor(2*N, -1);
    // 探索パートっぽいところ
    vector<int> done(N);
    for (int i = delta-1; ; i--) {
        for (int i = 0; i < N; i++)
            done[i] = 0;
        bool finish = true;
        vector<vector<int> > nD;
        double dist = beta * pow(2, i-1) * mini;
        double cost = pow(2, i+1) * mini;
        //cout << dist << " " << cost << endl;
        for (int l = 0; l < N; l++) {
            int num = D.size();
            for (int j = 0; j < num; j++) {
                if (D[j].size() == 1) {
                    cor[node-num+j] = D[j][0];
                    continue;
                }
                finish = false;
                vector<int> vs;
                for (int v : D[j]) {
                    if (d[perm[l]][v] <= dist && !done[v]) {
                        done[v] = 1;
                        vs.push_back(v);
                    }
                }
                if (vs.size()) {
                    T[node-num+j].emplace_back(node+nD.size(), cost);
                    nD.push_back(vs);
                }
            }
        }
        node += nD.size();
        D = nD;
        if (finish) break;
    }
    // 結果を出力
    for (int i = 0; i < node; i++) {
        cout << i << " " << cor[i] << endl;
        for (Edge e : T[i]) {
            cout << e.to << " " << e.cost << endl;
        }
        cout << endl;
    }
}