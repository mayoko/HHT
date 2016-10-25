// 木の情報を受け取って, 任意の 2 点間の距離を求めるプログラム
// O(N^2)

#include<iostream>
#include<vector>
#include<cstring>
#include<algorithm>

using namespace std;

struct Edge {
    int to;
    int cost;
    Edge() {}
    Edge(int t, int c) : to(t), cost(c) {}
};

// 頂点の最大数
const int MAXV = 5555;
// 木グラフ
vector<Edge> T[2*MAXV];
// 対応する頂点
int cor[2*MAXV];
int rcor[2*MAXV];
// 頂点の親
int par[2*MAXV];
int lca[2*MAXV][2*MAXV];
// 深さ(辺のコストは考慮してない)
int depth[2*MAXV];
// 辺のコストを考慮した深さ
int length[2*MAXV];
// 木のどこまでが実際の頂点と対応してるか的な
int maxID;

int ans[MAXV][MAXV];

void dfs(int v, int d, int l) {
    depth[v] = d;
    length[v] = l;
    for (Edge e : T[v]) {
        dfs(e.to, d+1, l+e.cost);
    }
}

int main() {
    // 入力
    int N;
    cin >> N;
    memset(par, -1, sizeof(par));
    for (int i = 0; i < N; i++) {
        int id;
        cin >> id >> cor[i];
        maxID = max(maxID, cor[i]);
        if (cor[i] != -1) rcor[cor[i]] = i;
        int edgeNum;
        cin >> edgeNum;
        T[i].resize(edgeNum);
        for (int j = 0; j < edgeNum; j++) {
            cin >> T[i][j].to >> T[i][j].cost;
            par[T[i][j].to] = i;
        }
    }
    // 深さを求める
    dfs(0, 0, 0);
    // lca を求める
    {
        vector<int> order(N);
        for (int i = 0; i < N; i++)
            order[i] = i;
        sort(order.begin(), order.end(), [&](int i, int j) {return depth[i] < depth[j];});
        for (int i = 0; i < N; i++) {
            int v = order[i];
            for (int j = 0; j < N; j++) {
                int u = order[j];
                if (v == u) lca[v][u] = v;
                else if (depth[v] > depth[u]) lca[v][u] = lca[par[v]][u];
                else if (depth[v] < depth[u]) lca[v][u] = lca[v][par[u]];
                else if (depth[v] > 0) lca[v][u] = lca[par[v]][u];
                else lca[v][u] = 0;
            }
        }
    }
    // 答えの出力
    for (int i = 0; i <= maxID; i++) {
        for (int j = 0; j <= maxID; j++) {
            int v = rcor[i], u = rcor[j];
            int l = lca[v][u];
            ans[i][j] = length[v]-length[l] + length[u]-length[l];
        }
    }
    for (int i = 0; i <= maxID; i++) {
        for (int j = 0; j <= maxID; j++) {
            cout << ans[i][j] << " " ;
        }
        cout << endl;
    }
}