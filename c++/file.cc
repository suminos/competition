// 標準入出力をファイル入出力に切り替える
int main() {
  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);
}

// fopen, fclose を使う
int main() {
  FILE *in = fopen("input.txt", "r");
  int n; fscanf(in, "%d", &n);
  fclose(in);
  
  FILE *out = fopen("output.txt", "w");
  fprintf(out, "%d\n", n);
  fclose(out);
}

