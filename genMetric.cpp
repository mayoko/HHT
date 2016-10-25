// 三角不等式を満たす距離空間を作るプログラム
// グリッド上に適当に点をばらまいてそれのマンハッタン距離を求める
// O(N^2 log N)
#include<iostream>
#include<random>
#include<vector>
#include<set>
#include<utility>

using namespace std;

random_device rnd;
mt19937 mt(rnd());

int main() {
    // グラフの頂点数
    int N;
    cin >> N;
    // 距離の最大値の目安(適当)
    int MAX;
    cin >> MAX;
    set<pair<int, int> > S;
    while (S.size() < N) {
        S.insert(make_pair(mt()%MAX, mt()%MAX));
    }
    vector<pair<int, int> > vs(S.begin(), S.end());
    // 返す距離空間は行列で作る
    vector<vector<int> > mat(N, vector<int>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            mat[i][j] = abs(vs[i].first-vs[j].first) + abs(vs[i].second-vs[j].second);
        }
    }
    // 出力
    cout << N << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }

    vector<vector<int> > req(N, vector<int>(N));
    // 後 requirement 用に適当に作っておく
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            req[i][j] = req[j][i] = mt()%MAX;
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << req[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}
