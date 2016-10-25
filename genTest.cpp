#include<bits/stdc++.h>
using namespace std;
random_device rnd;
mt19937 mt(rnd());


int main() {
    int N, MAX;
    cin >> N >> MAX;
    vector<vector<int> > req(N, vector<int>(N));
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            req[i][j] = req[j][i] = mt()%MAX;
        }
    }
    cout << N << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << req[i][j] << " ";
        }
        cout << endl;
    }
    cout << 10 << endl;
    return 0;
}
