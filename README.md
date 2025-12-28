# DSP Assignment 4: Sampling Rate Conversion with FFT Filters

##  內容摘要
本專案實作了數位訊號處理中的取樣率轉換 (Sampling Rate Conversion)。
目標將 **44.1kHz** 的立體聲音樂檔案 (`input.wav`) 轉換為 **8kHz** 的輸出檔案 (`output.wav`)。

## 系統規格
* **Input Signal**: 44.1 kHz, 16-bit PCM, Stereo (`input.wav`)
* **Output Signal**: 8.0 kHz, 16-bit PCM, Stereo (`output.wav`)
* **Language**: C 
* **Algorithm**: FFT-based FIR Filtering using Overlap-Add Method

## 系統參數
為了達成 $44100 \text{ Hz} \to 8000 \text{ Hz}$ 的轉換，我們計算出整數的升頻與降頻因子：

$$\text{Target Rate} = \text{Input Rate} \times \frac{L}{M}$$
$$8000 = 44100 \times \frac{80}{441}$$

因此系統參數設定如下：
* **Up-sampling factor ($L$)**: 80
* **Down-sampling factor ($M$)**: 441
* **Frame Size ($P$)**: 441 (Length of segment for processing)
* **Filter Length ($Q$)**: 1025
* **FFT Size ($N$)**: 2048 (Next power of 2 such that $N \ge P + Q - 1$)

## 編譯方式
由於程式包含大量的浮點數運算與迴圈，經過實測後如果使用正常`-lm` 參數的話要跑超過10分鐘，所以**強烈建議** 使用 `-O3` 參數進行編譯優化，以大幅提升執行速度。

### 編譯
```bash
gcc main.c -o main.exe -lm -O3
```

### 執行
```bash
main.exe
```
------

## 程式碼詳細解析
我會把程式碼拆解成七個部分，並解釋每一行背後的數學意義與程式邏輯。
### 第一部分：標頭檔與常數定義
``` C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
```
* `<math.h>`：為了使用 `sin`,` cos` (計算 FFT 的旋轉因子與濾波器形狀)。
* `PI`：定義圓周率，因為傅立葉變換與 Sinc 函數都需要用到 $\pi$。
* **投影片對應：** 這是為了實現 FFT 與濾波器公式的基礎建設。

----
### 第二部分：資料結構 
這部分定義了兩個結構，這是 C 語言處理「複數」與「WAV 檔案」的基礎。
``` C
typedef struct {
    double r; // Real part (實部)
    double i; // Imaginary part (虛部)
} Complex;
```
* **解釋：** 標準 C 語言沒有內建複數型態（除非用 C99 的 complex.h，但作業通常希望從底層寫起）。FFT 運算必須在複數域進行，所以我們定義 r 代表實部，i 代表虛部。
``` C
typedef struct {
    char riff[4];        // 檔案標記 "RIFF"
    unsigned int size;   // 檔案大小
    // ... (中間略過) ...
    unsigned int sample_rate; // 取樣率 (如 44100)
    // ...
    unsigned int data_size;   // 音訊數據的總位元組數
} WavHeader;
```
* 　**用途**：
    1. 讀檔時：我們讀取這個結構，知道音樂有多長、取樣率是多少。
    2. 寫檔時：我們必須修改 sample_rate 為 8000，並更新 data_size，否則播放器會無法播放。
----
### 第三部分：複數運算與 FFT 核心
#### 1.複數基礎運算
``` C
Complex mul(Complex a, Complex b) {
    Complex c;
    c.r = a.r * b.r - a.i * b.i; // 公式：(ac - bd)
    c.i = a.r * b.i + a.i * b.r; // 公式：i(ad + bc)
    return c;
}
```
* **解釋** ：這是複數乘法公式 $(a+bi)(c+di)$ 的實作。這在 頻域濾波 (卷積) 時會被大量使用。
#### 2.FFT 函式
```C
void fft(Complex *x, int n, int invert) { ... }
```
* `n`：FFT 的點數，這裡是 2048 。
* `invert`：`0` 代表做 FFT，`1` 代表做 IFFT (反轉換)。
#### 步驟A : 位元反轉 (Bit-reversal Permutation)
```C
int j = 0;
    for(int i = 0; i < n; i++) {
        if(i < j) {
            // 交換 x[i] 和 x[j]
            // ...
        }
        // ... (計算反轉索引 j 的邏輯)
    }
```
* **原理：** Cooley-Tukey 演算法需要將輸入陣列「打亂」，把偶數索引和奇數索引分開。經過位元反轉排序後，我們才能用迴圈由下而上進行合併運算。
#### 步驟B : 蝴蝶運算 (Butterfly Operation)
```C
for (int len = 2; len <= n; len <<= 1) { // 第一層：級數 (2, 4, 8, ..., 2048)
        double ang = 2 * PI / len * (invert ? -1 : 1); // 計算旋轉角度
        Complex wlen = {cos(ang), sin(ang)}; // 單位根 (旋轉因子)
        
        for(int i = 0; i < n; i += len) { // 第二層：遍歷每一個區塊
            Complex w = {1.0, 0.0};
            for(int k = 0; k < len / 2; k++) { // 第三層：蝴蝶翅膀
                // 核心運算：上 = u + v*w, 下 = u - v*w
                // ...
            }
        }
    }
```
* **解釋：** 這是將小點數的 DFT 組合成大點數 DFT 的過程。w (Twiddle Factor) 代表旋轉因子 $W_N^k$，是頻率分析的核心。
#### 步驟C : IFFT 正規化
```C
if(invert) {
        for(int i = 0; i < n; i++) {
            x[i].r /= n; // 除以 N
            // ...
        }
    }
```
* **解釋**：依據數學定義，做完反傅立葉變換後，數值會放大 $N$ 倍，必須除回去才能還原正確的振幅。

----

### 第四部分：濾波器設計 (Filter Design)
這部分實作了 低通濾波器，目的是消除升頻產生的鏡像雜訊。
```C
void design_lowpass(Complex *h_buffer, int N, int Q, int L, int M) {
    double cutoff = PI / M; // 截止頻率設定為 pi/441 [cite: 20]
    double center = (Q - 1) / 2.0; // 濾波器中心點，確保相位線性

    for (int n = 0; n < N; n++) {
        if (n < Q) {
            // 1. Sinc 函數：理想低通濾波器的時域形狀
            double val = sin(cutoff * x) / (PI * x);

            // 2. Blackman Window：修飾 Sinc 函數
            // 如果直接截斷 Sinc，頻譜會有強烈震盪 (Gibbs Phenomenon)。
            // 乘上 Blackman 窗函數可以讓邊緣平滑，大幅減少雜訊。
            double w = 0.42 - 0.5 * cos(...) + 0.08 * cos(...);
            val *= w;

            // 3. 增益補償 (Gain Compensation)
            val *= L; // 因為升頻補了很多 0，能量被稀釋了，需乘上 80 倍補回來 [cite: 14]。
            
            h_buffer[n].r = val;
        } else {
            // 補零 (Zero Padding)：因為 FFT 需要 2048 點，濾波器只有 1025 點 。
            h_buffer[n].r = 0.0; 
        }
    }
}
```
----

### 第五部分：主程式初始化
```C
const int L = 80;   // 升頻因子
const int M = 441;  // 降頻因子
const int P = 441;  // 每次處理的區塊長度 [cite: 44]
const int Q = 1025; // 濾波器長度 [cite: 45]
const int N = 2048; // FFT 點數 (必須 >= P + Q - 1 以避免混疊)
```
* **記憶體配置：**
    * `buffer`: 用來放音訊做 FFT 的暫存區。
    * `H`: 用來放濾波器的頻譜（算好一次就可以一直用）。
    * `overlap_buf`: 用來做 Overlap-Add 的累加緩衝區。
* **濾波器準備：**
    ```C
    design_lowpass(buffer, ...); // 產生時域濾波器 h[n]
    fft(buffer, N, 0);           // 轉成頻域 H[k]
    // 存入 H 陣列供後續使用
    ```
    
---

### 第六部分：核心迴圈
這段程式碼模擬了訊號經過系統的過程。
```C
while (current_pos < total_xe_len){
```
*` current_pos`：這是一個「虛擬」的時間軸，代表升頻後的訊號位置。
#### 步驟A : 虛擬升頻 (Virtual Up-sampling)
```C
for(int k = 0; k < N ; k++){
            // ...
            if(k < P) {
                long long xe_idx = current_pos + k;
                // 判斷是否為原始樣本點 (每 L 點出現一次)
                if(xe_idx % L == 0) { 
                    int original_idx = (int) (xe_idx / L);
                    // 填入原始資料，其他位置保持 0 (補零)
                    buffer[k].r = (double)input_data[original_idx];
                }
            }
        }
```
* **重點：** 我們沒有真的建立一個巨大的陣列來放升頻後的資料（那樣會塞爆記憶體）。我們用` if (xe_idx % L == 0)` 動態判斷：如果是 80 的倍數就填入資料，否則就是 0。這完全符合 **Up-sampling** 的數學定義 。

#### 步驟 B: 快速卷積 (FFT Convolution)
```C
fft(buffer, N, 0); // 1. 時域 -> 頻域
for(...) buffer[i] = mul(buffer[i], H[i]); // 2. 頻域相乘 (相當於時域卷積) 
fft(buffer, N, 1); // 3. 頻域 -> 時域
```
* 這是 DSP 最標準的加速技巧。直接做卷積太慢，轉去頻域用乘的比較快。
#### 步驟 C: 重疊相加 (Overlap-Add)
```C
for(int i = 0; i < N; i++) {
    overlap_buf[i] += buffer[i].r;
}
```
* 因為濾波器有長度，每個區塊運算完會有「尾巴」(Transient)。這個尾巴必須留給下一個區塊加在一起，聲音才會連續。這就是投影片要求的 **Overlap-Add method** 。
#### 步驟 D: 降頻與輸出 (Down-sampling)
```
if(out_idx < output_capacity) {
    double val = overlap_buf[0]; // 只取第一個點！
    output_data[out_idx++] = (short)(val);
}
```
* 為什麼只取 overlap_buf[0]？
    * 我們每次處理 $P=441$ 個點。 
    * 我們要降頻 $M=441$，也就是每 441 個點留 1 個。
    * 這是一個精巧的設計：每一輪迴圈結束，剛好產出一個我們要的點。其他的點雖然算出來了，但 在降頻過程中被「丟棄」了。
#### 步驟 E: 位移緩衝區 (Shift)
```C
for (int i = 0; i < N - P; i++) {
    overlap_buf[i] = overlap_buf[i + P];
}
```
* 把舊的資料推出去，把位置讓出來給下一輪的運算結果。

---

#### 第七部分：寫檔與收尾
```C
header.sample_rate = 8000; // 修改 Header 告訴播放器這是 8kHz
// ...
fwrite(...); // 寫入硬碟
free(...);   // 釋放記憶體，好習慣
```
---
## 流程總結
整個程式的運作流程如下：
1. 讀檔 (`input.wav`)
2. 準備濾波器 ($H$)：計算 Sinc + Window $\rightarrow$ 做一次 FFT 存起來。
3. 迴圈開始：
    * 抓取資料：從 `input` 抓一點點資料，中間補 79 個 0 (邏輯上)。
    * FFT：變成頻譜。
    * 乘法：頻譜 $\times H$ (把高頻雜訊切掉)。
    * IFFT：變回波形。
    * 累加：把波形加到 `overlap_buffer`。
    * 輸出：從 `overlap_buffer` 拿第 0 個點存檔。
    * 位移：把` overlap_buffer` 往前推。
4. 寫檔 (`output.wav`)