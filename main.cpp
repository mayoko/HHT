#include<iostream>
#include<random>
#include<vector>
#include<algorithm>
#include<unordered_map>
#include<cassert>
#include<cstring>
#include<queue>

using namespace std;

random_device rnd;
mt19937 mt(rnd());
// 距離の最大値
const int MAX = 114514;
// 頂点数の最大値
const int MAXV = 3333;
// 頂点数(入力)
int N;

struct Edge {
    int to;
    int cost;
    int cap;
    Edge() {}
    Edge(int t, int c) : to(t), cost(c) {}
};

// 距離空間
int dist[MAXV][MAXV];
// requirement
int R[MAXV][MAXV];

// 木(出力すべきもの)
vector<Edge> T[2*MAXV];
// 木の頂点に対応する頂点群がシングルトンになっているとき, 対応する頂点をメモしておく(これも出力すべきもの)
int corr[2*MAXV];
// 今のところの木の頂点数
int treeNum;

// 距離空間が与えられるので, 近似的に木距離にする
// O(N^2)
namespace approxTreeMetric {
    // 木が何段組みになるかってヤツ
    int delta;
    // 例の乱数
    double beta;
    // 例のランダム順列
    int perm[MAXV];
    // 距離空間で最大の距離と最小の距離
    int mini = MAX, maxi = 0;
    // ある頂点が, 今順列の何番目まで進んでいるか
    int cur[MAXV];

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
        double radius = beta * pow(2, delta-depth-2)*mini;
        int cost = pow(2, delta-depth)*mini;
        for (int el : vs) {
            int& now = cur[el];
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

    // 与えられた入力を用いて解く
    void solve() {
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
        cerr << treeNum << endl;
        for (int i = 0; i < treeNum; i++) {
            cerr << i << " " << corr[i] << endl;
            for (Edge e : T[i]) {
                cerr << e.to << " " << e.cost << endl;
            }
            cerr << endl;
        }
    }
}

// 木のどこまでが実際の頂点と対応してるか的な
int maxID;
// 木グラフから距離に変換する
namespace treeToMetric {
    // 出力すべき距離
    int ans[MAXV][MAXV];
    // 対応する頂点
    int rcorr[2*MAXV];
    // 頂点の親
    int par[2*MAXV];
    int lca[2*MAXV][2*MAXV];
    // 深さ(辺のコストは考慮してない)
    int depth[2*MAXV];
    // 辺のコストを考慮した深さ
    int length[2*MAXV];

    void dfs(int v, int d, int l) {
        depth[v] = d;
        length[v] = l;
        for (Edge e : T[v]) {
            dfs(e.to, d+1, l+e.cost);
        }
    }

    void solve() {
        int N = treeNum;
        for (int i = 0; i < N; i++) {
            maxID = max(maxID, corr[i]);
            if (corr[i] != -1) rcorr[corr[i]] = i;
            int edgeNum = T[i].size();
            for (int j = 0; j < edgeNum; j++) {
                par[T[i][j].to] = i;
            }
        }
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
        for (int i = 0; i <= maxID; i++) {
            for (int j = 0; j <= maxID; j++) {
                int v = rcorr[i], u = rcorr[j];
                int l = lca[v][u];
                ans[i][j] = length[v]-length[l] + length[u]-length[l];
                //cerr << ans[i][j] << " " ;
            }
            //cerr << endl;
        }
        vector<double> muls;
        double avg = 0, mini = MAX, maxi = 0, sum = 0, total;
        for (int i = 0; i <= maxID; i++) {
            for (int j = i+1; j <= maxID; j++) {
                total += ans[i][j];
                double tmp = 1. * ans[i][j] / dist[i][j];
                assert(tmp >= 1);
                mini = min(mini, tmp);
                maxi = max(maxi, tmp);
                muls.push_back(tmp);
                avg += tmp;
                sum += dist[i][j];
            }
        }
        sort(muls.begin(), muls.end());
        cerr << "total ratio: " << total/sum << endl;
        int num = (maxID+1)*maxID/2;
        avg /= num;
        cerr << "avg ratio: " << avg << endl;
        cerr << "min ratio: " << mini << endl;
        cerr << "med ratio: " << muls[num/2] << endl;
        cerr << "max ratio: " << maxi << endl;
    }
}

// 与えられた木を trivalent tree に変換する
// trivalent tree は無向グラフにしておく
// ところでこのアルゴリズムは葉の頂点数が 3 以上であることを前提にしていることに注意
// O(N)

// 出力すべきもの
vector<Edge> Tri[2*MAXV];
namespace makeTrivalentTree {
    // 一時的に作る木
    vector<Edge> tT[2*MAXV];
    // 木からこの頂点を除去するかどうか
    bool remove[2*MAXV];
    // 各頂点の次数
    int degree[2*MAXV];
    // 各頂点の深さ
    int depth[2*MAXV];

    inline bool isValid(const Edge& e, const int p) {
        return e.to != p && !remove[e.to];
    }
    void getDepth(int v, int p, int d) {
        depth[v] = d;
        for (Edge e : tT[v]) if (isValid(e, p)) {
            getDepth(e.to, v, d+e.cost);
        }
    }

    // 今いる頂点は v
    // 親は p
    // パスで続いている時その一端は pre
    void dfs(int v, int p, int pre) {
        if (degree[v] == 1) {
            // 葉なので pre とつなげる
            Tri[pre].emplace_back(v, depth[v] - depth[pre]);
            return;
        } else if (degree[v] == 2) {
            // ただのパスなのでスルー
            for (Edge e : tT[v]) if (isValid(e, p)) {
                dfs(e.to, v, pre);
            }
            return;
        }
        // pre からの系譜は終了したので辺を張る
        if (pre != -1) Tri[pre].emplace_back(v, depth[v]-depth[pre]);
        if (degree[v] == 3) {
            // Tri に辺を追加
            for (Edge e : tT[v]) if (isValid(e, p)) {
                dfs(e.to, v, v);
            }
        } else {
            int addNum = 1;
            // root の場合は頂点に直属の頂点が増える
            if (p == -1) ++addNum;
            // とりあえず連結リストっぽいやつを作る
            Tri[v].emplace_back(treeNum, 0);
            depth[treeNum] = depth[v];
            degree[treeNum] = 3;
            // 新しく追加する頂点の数は degree[v]-3 個
            for (int i = 0; i < degree[v]-3-1; i++) {
                Tri[treeNum+i].emplace_back(treeNum+i+1, 0);
                depth[treeNum+i+1] = depth[v];
                degree[treeNum+i+1] = 3;
            }
            // v の子リスト
            vector<int> list;
            for (int i = 0; i < tT[v].size(); ++i) if (isValid(tT[v][i], p)) {
                list.push_back(i);
            }
            assert(list.size() == degree[v]-(p!=-1));
            int tn = treeNum;
            treeNum += degree[v]-3;
            // v 直属おじさん
            for (int i = 0; i < addNum; ++i) {
                Edge e = tT[v][list[i]];
                dfs(e.to, v, v);
            }
            // それ以外
            for (int i = addNum; i < list.size(); ++i) {
                Edge e = tT[v][list[i]];
                int num = tn+i-addNum;
                if (i == list.size()-1) --num;
                dfs(e.to, v, num);
            }
            //for (int i = 0; i < addNum; ++i) {
            //    Edge e = tT[v][list[i]];
            //    Tri[v].emplace_back(e.to, e.cost);
            //}
            //for (int i = addNum; i < list.size(); ++i) {
            //    Edge e = tT[v][list[i]];
            //    Tri[treeNum+i-addNum].emplace_back(e.to, e.cost);
            //}
        }
    }

    void solve() {
        // approxTreeMetric で作った木は有向木になっているので無向木にする
        for (int i = 0; i < treeNum; i++) {
            for (Edge e : T[i]) {
                tT[i].emplace_back(e.to, e.cost);
                tT[e.to].emplace_back(i, e.cost);
            }
        }
        // 次数が 1 のものでもとの頂点でないやつは関係ないので消す
        // トポロジカルソートの要領
        queue<int> leaves;
        for (int i = 0; i < treeNum; i++) {
            degree[i] = tT[i].size();
            if (corr[i] < 0 && degree[i] == 1) {
                remove[i] = true;
                leaves.push(i);
            }
        }
        while (!leaves.empty()) {
            int now = leaves.front(); leaves.pop();
            for (Edge e : tT[now]) {
                int ch = e.to;
                if (degree[ch] == 1) continue;
                if (--degree[ch] == 0) {
                    // 元の頂点は消されないはず
                    assert(corr[ch] < 0);
                    remove[ch] = true;
                    leaves.push(ch);
                }
            }
        }
        // 次数が 3 以上の頂点を探す
        int start = -1;
        for (int i = 0; i < treeNum; i++) {
            if (degree[i] >= 3) {
                start = i;
                break;
            }
        }
        // 頂点数が 3 以上なら必ず start が存在する
        assert(start != -1);
        getDepth(start, -1, 0);
        // 木 dfs して trivalent tree を構成する
        dfs(start, -1, -1);

        // 余計な頂点があるので消しておく
        vector<int> A;
        // 必要な頂点を A にいれていく
        for (int i = 0; i < treeNum; i++) {
            tT[i] = Tri[i];
            if (corr[i] >= 0 || Tri[i].size() == 2 || Tri[i].size() == 3) {
                A.push_back(i);
            }
            Tri[i].clear();
        }
        for (int i = 0; i < A.size(); ++i) {
            // corr の整理
            if (corr[A[i]] >= 0) {
                corr[i] = corr[A[i]];
            }
            // Tri の整理
            for (Edge e : tT[A[i]]) {
                int num = lower_bound(A.begin(), A.end(), e.to)-A.begin();
                Tri[i].emplace_back(num, e.cost);
                Tri[num].emplace_back(i, e.cost);
            }
        }
        // もう T, tT はいらないので clear しておく
        for (int i = 0; i < treeNum; i++) {
            tT[i].clear();
            T[i].clear();
        }
        treeNum = A.size();
        for (int i = 0; i < treeNum; i++) {
            if (Tri[i].size() == 1) {
                assert(corr[i] >= 0);
            } else {
                assert(Tri[i].size() == 3);
            }
        }

        // デバッグ: 距離があっているか
        //const int INF = 1e9;
        //vector<vector<int> > d(treeNum, vector<int>(treeNum, INF));
        //for (int i = 0; i < treeNum; ++i) if (!remove[i]) {
        //    for (Edge e : Tri[i]) {
        //        d[i][e.to] = d[e.to][i] = e.cost;
        //    }
        //    d[i][i] = 0;
        //}
        //for (int k = 0; k < treeNum; k++) {
        //    for (int i = 0; i < treeNum; i++) {
        //        for (int j = 0; j < treeNum; j++) {
        //            d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
        //        }
        //    }
        //}
        //vector<vector<int> > hope(maxID+1, vector<int>(maxID+1));
        //for (int i = 0; i < treeNum; i++) if (corr[i] >= 0) {
        //    for (int j = 0; j < treeNum; j++) if (corr[j] >= 0) {
        //        int v = corr[i], u = corr[j];
        //        hope[v][u] = d[i][j];
        //    }
        //}
        //for (int i = 0; i <= maxID; i++) {
        //    for (int j = 0; j <= maxID; j++) {
        //        assert(hope[i][j] == treeToMetric::ans[i][j]);
        //    }
        //}

        // デバッグ出力
        cerr << treeNum << endl;
        for (int i = 0; i < treeNum; i++) if (!remove[i]) {
            cerr << i << " " << corr[i] << endl;
            for (Edge e : Tri[i]) {
                cerr << e.to << " " << e.cost << endl;
            }
            cerr << endl;
        }
    }
}

// N 頂点間の requirement R および木 T が与えられるので, 以下のクエリに答える
// 木 T のある辺で木を集合 X, ,Y の二つに分割するので, i \in X, j \in Y であるようなもので最大の R[i][j] を求める
// このクエリに答えるには, 実は R の最大全域木を知っておけば十分
// 集合 X のクエリに答えるには, 以下のようにすれば良い:
// 求めた最大全域木について, X 側に 1 つだけ端点を持つ辺の requirement の最大値を選ぶ
// これは O(N) なので神
namespace Oracle {
    // 最大全域木
    vector<int> dominantTree[MAXV];
    // assign されたかどうか
    bool done[MAXV];
    // assign されてない頂点について, その親の候補
    int par[MAXV];
    void makeDominantTree() {
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
            dominantTree[par[v]].push_back(v);
            dominantTree[v].push_back(par[v]);
            done[v] = true;
            for (int i = 0; i < N; i++) if (!done[i]) {
                if (R[par[i]][i] < R[v][i]) par[i] = v;
            }
        }
        for (int i = 0; i < N; i++) {
            cerr << i << endl;
            for (int el : dominantTree[i]) cerr << el << " ";
            cerr << endl;
        }
    }
    bool isIn[MAXV];
    // 頂点集合に対して, 解を出力する
    int calc(const vector<int>& vs) {
        memset(isIn, false, sizeof(isIn));
        for (int el : vs) isIn[el] = true;
        int ret = 0;
        for (int v : vs) {
            for (int u : dominantTree[v]) {
                if (!isIn[u]) ret = max(ret, R[v][u]);
            }
        }
        return ret;
    }
}

int solutionValue = 0;
// 作った trivalent tree Tri に対して, Edge::cap を求めていく
namespace calcEdgeWeight {
    vector<int> dfs(int v, int p) {
        vector<int> ret;
        bool leaf = true;
        for (Edge& e : Tri[v]) if (e.to != p) {
            leaf = false;
            vector<int> tmp = dfs(e.to, v);
            e.cap = Oracle::calc(tmp);
            for (Edge& f : Tri[e.to]) {
                if (f.to == v) f.cap = e.cap;
            }
            ret.insert(ret.end(), tmp.begin(), tmp.end());
        }
        if (leaf) {
            assert(corr[v] >= 0);
            ret.push_back(corr[v]);
        }
        cerr << v << endl;
        for (int el : ret) cerr << el << " ";
        cerr << endl << endl;
        return ret;
    }
    void solve() {
        cerr << "start dfs" << endl;
        dfs(0, -1);
        // デバッグ出力
        cerr << treeNum << endl;
        for (int i = 0; i < treeNum; i++) {
            cerr << i << " " << corr[i] << endl;
            for (Edge e : Tri[i]) {
                cerr << e.to << " " << e.cost << " " << e.cap << endl;
            }
            cerr << endl;
        }
    }
}

vector<Edge> ans[MAXV];
namespace getSolution {
    void setTri(int v, int p) {
        for (Edge& e : Tri[v]) {
            e.cap *= 2;
            if (e.to != p) setTri(e.to, v);
        }
    }

    int dfs(int v, int p, const Edge& e, int s, int mini) {
        cerr << v << endl;
        if (Tri[v].size() == 1) {
            assert(mini != 0);
            ans[corr[s]].emplace_back(corr[v], (mini+1)/2);
            ans[corr[v]].emplace_back(corr[s], (mini+1)/2);
            for (Edge& te : Tri[v]) {
                if (te.to == p) te.cap -= mini;
            }
            return mini;
        }
        vector<int> valid(2);
        int cnt = 0, pIndex;
        for (int i = 0; i < 3; i++) {
            Edge te = Tri[v][i];
            if (te.to != p) valid[cnt++] = i;
            else pIndex = i;
        }
        int index = !(e.cap + Tri[v][valid[0]].cap - Tri[v][valid[1]].cap > 0);
        assert(e.cap + Tri[v][valid[index]].cap - Tri[v][valid[index^1]].cap > 0);
        Edge e0 = Tri[v][valid[index]];
        int tmp = (e.cap + e0.cap - Tri[v][valid[index^1]].cap)/2;
        mini = min(mini, tmp);
        int ret = dfs(e0.to, v, e0, s, mini);
        Tri[v][valid[index]].cap -= ret;
        Tri[v][pIndex].cap -= ret;
        return ret;
    }

    void solve() {
        // Tri 木の cap をすべて 2 倍する
        setTri(0, -1);
        // すべての葉の cap が 0 になるまでパスを作りつづける
        const int INF = 1e9;
        while (1) {
            bool finish = true;
            for (int i = 0; i < treeNum; ++i) {
                if (Tri[i].size() == 1) {
                    assert(corr[i] >= 0);
                    Edge& e = Tri[i][0];
                    if (e.cap > 0) {
                        finish = false;
                        cerr << "start " << i << endl;
                        int tmp = dfs(e.to, i, e, i, INF);
                        cerr << endl;
                        e.cap -= tmp;
                    }
                }
            }
            if (finish) break;
        }
        // 解を出力
        for (int i = 0; i < N; ++i) {
            for (Edge e : ans[i]) {
                solutionValue += e.cost * dist[i][e.to];
                cout << i << " " << e.to << " " << e.cost << endl;
            }
        }
        cerr << "solution value is: " << solutionValue << endl;
    }
}

int main() {
    cin >> N;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> dist[i][j];
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> R[i][j];
        }
    }
    // 近似する
    cerr << "############################" << endl;
    cerr << "approx Tree Metric" << endl;
    approxTreeMetric::solve();
    // 木距離を求める
    cerr << "############################" << endl;
    cerr << "tree to metric" << endl;
    treeToMetric::solve();
    // trivalent tree に変換
    cerr << "############################" << endl;
    cerr << "make trivalent tree" << endl;
    makeTrivalentTree::solve();
    // oracle のための最大全域木を作る
    cerr << "############################" << endl;
    cerr << "make dominant tree" << endl;
    Oracle::makeDominantTree();
    // edgeWeight を求めて Tri に入れていく
    cerr << "############################" << endl;
    cerr << "calc Edge Weight" << endl;
    calcEdgeWeight::solve();
    // 感動の最終話
    cerr << "############################" << endl;
    cerr << "get solution" << endl;
    getSolution::solve();
}
