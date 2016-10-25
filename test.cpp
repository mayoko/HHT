#include<bits/stdc++.h>

using namespace std;
random_device rnd;
mt19937 mt(rnd());
const int MAXN = 555;
int R[MAXN][MAXN];
int N;

bool isIn[MAXN];
int naive(vector<int> vs) {
    memset(isIn, false, sizeof(isIn));
    for (int el : vs) isIn[el] = true;
    int ret = 0;
    for (int i = 0; i < N; i++) if (isIn[i]) {
        for (int j = 0; j < N; j++) if (!isIn[j]) {
            ret = max(ret, R[i][j]);
        }
    }
    return ret;
}

int solve(const vector<int>& vs, const vector<vector<int> >& G) {
    memset(isIn, false, sizeof(isIn));
    for (int el : vs) isIn[el] = true;
    int ret = 0;
    for (int v : vs) {
        for (int u : G[v]) {
            if (!isIn[u]) ret = max(ret, R[v][u]);
        }
    }
    return ret;
}

bool done[MAXN];
int par[MAXN];
vector<vector<int> > makeDominantTree() {
    vector<vector<int> > G(N);
    done[0] = true;
    for (int i = 0; i < N; i++) {
        par[i] = 0;
    }
    for (int t = 0; t < N-1; t++) {
        int v = -1;
        for (int i = 0; i < N; i++) if (!done[i]) {
            if (v == -1 || R[i][par[i]] > R[v][par[v]]) v = i;
        }
        assert(v != -1);
        G[par[v]].push_back(v);
        G[v].push_back(par[v]);
        done[v] = true;
        for (int i = 0; i < N; i++) if (!done[i]) {
            if (R[par[i]][i] < R[v][i]) par[i] = v;
        }
    }
    return G;
}

int main() {
    cin >> N;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> R[i][j];
        }
    }
    auto G = makeDominantTree();
    int Q;
    cin >> Q;
    while (Q--) {
        set<int> S;
        for (int t = 0; t < 2*N/3; t++) {
            S.insert(mt()%N);
        }
        vector<int> vs(S.begin(), S.end());
        cout << vs.size() << endl;
        if (solve(vs, G) == naive(vs)) cout << "SUCCESS!" << endl;
        else cout << "FAIL.." << endl;
    }
    return 0;
}
