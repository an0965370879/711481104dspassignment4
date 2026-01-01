#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

// ==========================================
// Part 1: 資料結構定義
// ==========================================

typedef struct {
    double r;
    double i;
} Complex;

typedef struct {
    char riff[4];
    unsigned int size;
    char wave[4];
    char fmt[4];
    unsigned int fmt_size;
    unsigned short format;
    unsigned short channels;
    unsigned int sample_rate;
    unsigned int byte_rate;
    unsigned short block_align;
    unsigned short bit_per_sample;
    char data[4];
    unsigned int data_size;
} WavHeader;

// ==========================================
// Part 2: 數學運算函式 (手刻 FFT)
// ==========================================

Complex add(Complex a, Complex b) {
    Complex c = {a.r+b.r, a.i+b.i};
    return c;
}
Complex sub(Complex a, Complex b) {
    Complex c = {a.r-b.r, a.i-b.i};
    return c;
}
Complex mul(Complex a, Complex b) {
    Complex c;
    c.r = a.r * b.r - a.i * b.i;    // 公式：(ac - bd)
    c.i = a.r * b.i + a.i * b.r;    // 公式：i(ad + bc)
    return c;
}   //複數乘法公式 (a+bi)(c+di) = (ac-bd) + i(ad+bc)

// 核心 FFT 函式
// x: 資料陣列, n: 點數, invert: 0=FFT, 1=IFFT

void fft(Complex *x, int n, int invert){
    
   int j = 0;
   for(int i = 0;i < n; i++) {
       if(i < j) {
           Complex temp = x[i];
           x[i] = x[j];     // 交換x[i]和x[j]
           x[j] = temp;     
       }
       int m = n / 2;
       while(m >= 1 && j >= m) { // 計算反轉索引j 的邏輯
           j -= m;
           m /= 2;
       }
       j += m;
   }


    // Butterfly Operations
    for (int len = 2; len<=n; len <<=1){ //// 第一層：級數 (2, 4, 8, ..., 2048)
        double ang = 2 * PI / len * (invert ? -1 : 1); // 計算旋轉因子角度
        Complex wlen = {cos(ang), sin(ang)}; // e^(i*ang) = cos(ang) + i*sin(ang)
        for(int i = 0; i < n; i += len){      // 第二層：每個子 DFT 的起始位置
            Complex w = {1.0, 0.0};
             for(int k = 0; k < len / 2; k++){   // 第三層：每個蝶形運算
                Complex u = x[i + k];
                Complex v = mul(x[i + k + len / 2], w); 
                x[i + k] = add(u, v);       // 上半部
                x[i + k + len / 2] = sub(u, v); // 下半部
              w = mul(w, wlen);
            }
        }
    }
    // IFFT Scaling
    if(invert) {
        for(int i = 0; i < n; i++) {
            x[i].r /= n;    // 對實部做縮放
            x[i].i /= n;    // 對虛部做縮放
        }
    }
}

// ==========================================
// Part 3: 設計濾波器
// ==========================================


// 產生低通濾波器並直接填入 FFT 輸入陣列
void design_lowpass(Complex *h_buffer, int N, int Q, int L, int M) {
    double cutoff = PI / M; // Cutoff = pi/441
    double center = (Q - 1) / 2.0;

    for (int n = 0; n < N; n++) {
        if (n < Q) {
            // 1. Sinc 函數：理想低通濾波器的時域形狀
            double val;
            double x = n - center;
            if (x == 0) val = cutoff / PI;
            else val = sin(cutoff * x) / (PI * x);

            // 2. Blackman Window：修飾 Sinc 函數
            // 如果直接截斷 Sinc，頻譜會有強烈震盪 (Gibbs Phenomenon)。
            // 乘上 Blackman 窗函數可以讓邊緣平滑，大幅減少雜訊。
            double w = 0.42 - 0.5 * cos(2 * PI * n / (Q - 1)) + 0.08 * cos(4 * PI * n / (Q - 1)); 
            val *= w;

            // 3. 增益補償 (Gain Compensation)
            val *= L; // 因為升頻補了很多 0，能量被稀釋了，需乘上 80 倍補回來 [cite: 14]。

            h_buffer[n].r = val; // 實部
            h_buffer[n].i = 0.0; // 虛部    
        } else {
            // 補零 (Zero Padding)：因為 FFT 需要 2048 點，濾波器只有 1025 點 。
            h_buffer[n].r = 0.0; // 實部
            h_buffer[n].i = 0.0; // 虛部
        }
    }
}


// ==========================================
// Part 4: 主程式
// ==========================================

int main(){

    int percent = 0; // 用來記錄百分比
    // 參數定義
    const int L = 80;
    const int M = 441;
    const int P = 441; // Processing block size (of upsampled signal)
    const int Q = 1025;// Filter length
    const int N = 2048;// FFT size (Power of 2, >= P+Q-1)

    printf("Starting fft \n");
    // 1. 讀取檔案
    FILE *fin = fopen("input.wav", "rb");
    if(!fin){
        printf("Error opening input.wav\n");
        return 1;
    }

    WavHeader header;
    fread(&header, sizeof(WavHeader), 1, fin);
    
    int num_input_samples = header.data_size / (header.bit_per_sample / 8);
    // 假設是單聲道或立體聲，這裡為了簡化示範，將資料視為一連串的 short
    // 若是立體聲，標準 WAV 是 L, R, L, R 交錯，此程式碼會依序處理 L 和 R (視為單一長串流)
    // 但正確的 DSP 應該 L 和 R 分開處理。
    // *注意*：若要嚴謹處理立體聲，建議把左右聲道拆開分別跑一次流程。
    // 以下程式碼將整個資料流視為單一訊號處理 (Interleaved processing works for linear filtering)

    short *input_data = (short *)malloc(num_input_samples * sizeof(short));
    fread(input_data, sizeof(short), num_input_samples, fin);
    fclose(fin);

    // 2. 準備記憶體
    Complex *H = (Complex *)malloc(N* sizeof(Complex)); // 濾波器頻譜
    Complex *buffer = (Complex *)malloc(N* sizeof(Complex)); // fft 工作區
    double *overlap_buf = (double *)calloc(N + P, sizeof(double)); // 重疊區域

    // 3. 製作並轉換濾波器
    design_lowpass(buffer, N, Q, L, M);
    fft(buffer, N, 0); // FFT
    for(int i = 0; i < N; i++)H[i] = buffer[i]; // 儲存濾波器頻譜

    // 4. 準備輸出
    // 輸出樣本數估算: Input * L / M
    int output_capacity = (int)((long long)num_input_samples*L / M) + 1000;
    short *output_data = (short *)malloc(output_capacity * sizeof(short));
    int out_idx = 0;
    // 5. 核心迴圈 (Overlap-Add)
    // 我們遍歷 "Up-sampled" 的虛擬時間軸
    long long total_xe_len = (long long)num_input_samples*L;
    long long current_pos = 0;

    while (current_pos < total_xe_len){
        // 每前進 1% 更新一次，避免頻繁 printf 拖慢速度
        if (current_pos * 100 / total_xe_len > percent) {
            percent = current_pos * 100 / total_xe_len;
            printf("Processing: %d%% \r", percent); // \r 會讓游標回到行首更新數字
        }
        // --- Step A: 準備 Segment (補零後的 Up-sample 訊號) ---
        for(int k = 0;k < N ;k++){
            buffer[k].r = 0.0;
            buffer[k].i = 0.0;

            // 只需要填前 P 個點
            if(k < P) {
                long long xe_idx = current_pos + k;
                // 檢查這點是否對應到原始輸入訊號(UP-sample邏輯)
                if(xe_idx % L == 0) {
                    int original_idx =(int) (xe_idx / L);
                    if(original_idx < num_input_samples) {
                        buffer[k].r = (double)input_data[original_idx];
                    }
                }
            }
        }
        // --- Step B: 卷積 (FFT -> Mult -> IFFT) ---

        fft(buffer, N, 0); // FFT
        for(int i = 0; i < N; i++) {
            buffer[i] = mul(buffer[i], H[i]); // 頻域相乘
        }
        fft(buffer, N, 1); // IFFT

        // --- Step C: Overlap-Add ---
        for(int i = 0; i < N; i++) {
            overlap_buf[i] += buffer[i].r; // IFFT 已經在函式內做過 /N 了
        }
        // --- Step D: Down-sampling (抽取) ---
        // 關鍵邏輯：
        // 我們的迴圈步長是 P = 441。
        // 我們的降頻因子是 M = 441。
        // 輸出 y[m] 對應時間點 m*M。
        // 當 m=0, 時間=0。當 m=1, 時間=441。
        // 剛好每次迴圈處理完，buffer 的第 0 個位置就是我們要的取樣點！

        if(out_idx < output_capacity) {
            // 取出 overlap_buf[0] 並轉回 short
            double val = overlap_buf[0];

            // Clipping (防爆音)
            if(val > 32767.0) val = 32767.0;
            if(val < -32768.0) val = -32768.0;
            output_data[out_idx++] = (short)(val);
        }
        // --- Step E: Shift Buffer ---
        // 移掉前 P 個點，後面的補上來
        // 這裡 P=441
        for (int i = 0; i < N - P; i++) {
            overlap_buf[i] = overlap_buf[i + P];
        }
        // 清除後面空出來的區域
        for (int i = N - P; i < N + P; i++) {
            overlap_buf[i] = 0.0;
        }
        current_pos += P; // 移動到下一個區塊

    }
    // 6. 寫入檔案
    FILE *fout = fopen("output.wav", "wb");

    // 更新 header 資訊
    header.sample_rate = 8000;
    header.byte_rate = 8000 * header.channels * (header.bit_per_sample / 8);
    header.data_size = out_idx * sizeof(short);
    header.size = header.data_size + 36;

    fwrite(&header, sizeof(WavHeader), 1, fout);
    fwrite(output_data, sizeof(short), out_idx, fout);
    fclose(fout);

    // 清理記憶體
    free(input_data);
    free(output_data);
    free(H);
    free(buffer);
    free(overlap_buf);

    printf("Finished fft: %d samples generated.\n", out_idx);
    return 0;
}
