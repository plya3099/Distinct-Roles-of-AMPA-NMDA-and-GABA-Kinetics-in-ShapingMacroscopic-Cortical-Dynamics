% 这是一个独立的MATLAB脚本，用于分析单个CSV文件。
% 功能：读取CSV数据，计算功率谱密度（PSD），排除DC分量，绘制频谱图，
% 计算Gamma频段（30-100 Hz）的归一化功率，以及主频（PSD最大值对应的频率）。
% 输出：保存PSD图像，并打印Gamma归一化功率和主频。
% 使用前，确保安装Signal Processing Toolbox（用于plomb和hann等函数）。

% 指定CSV文件路径（请修改为您的实际路径）
csv_path = 'D:\work\data\1114\14nmda\row_286\data.csv';%'D:\work\data\1103\row_8\data.csv';  % 示例路径，替换为您的文件

% 指定输出目录（可选，用于保存图像）
output_dir = fileparts(csv_path);  % 默认使用CSV所在目录
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 读取 CSV 文件
try
    df = readtable(csv_path);
catch
    error('无法读取 %s: 文件为空或读取错误', csv_path);
end

% 如果表格为空，停止
if height(df) == 0
    error('%s: 数据为空', csv_path);
end

% 假设列名为 't' 和 've'，如果不同，请修改
if ~ismember('t', df.Properties.VariableNames) || ~ismember('ve', df.Properties.VariableNames)
    error('%s: 缺少 ''t'' 或 ''ve'' 列', csv_path);
end

t = df.t * 0.02;  % 转换为秒，基于 dt=0.02 s（如果 df.t 已为秒，请移除 *0.02）
v = df.ve;

% 如果数据点少于 2，无法进行分析，停止
if length(t) < 2
    error('%s: 数据点不足', csv_path);
end

% 动态计算 dt
diffs = diff(t);
dt = mean(diffs);

% 检查数据均匀性（如果标准差过大，使用 plomb）
std_diff = std(diffs);
threshold = 0.005 * dt;  % 0.5% 相对阈值，可调整
use_plomb = std_diff > threshold;
if use_plomb
    fprintf('警告: 时间间隔不均匀 (std=%.6f > 阈值%.6f)，切换到Lomb-Scargle方法\n', std_diff, threshold);
end

% 调试：打印前几个 diff(t)
if length(t) >= 5
    fprintf('调试: 前几个 dt: %.8f, %.8f, %.8f, %.8f, %.8f\n', diffs(1:5));
else
    fprintf('调试: dt 差异: %s\n', mat2str(diffs, 8));
end

% 添加 NaN/Inf 检查
if any(isnan(v) | isinf(v))
    error('%s: 数据包含 NaN 或 Inf', csv_path);
end

% 去趋势（移除 DC 偏移）
v = detrend(v);

N = length(v);

if use_plomb
    % 使用 plomb 处理不均匀采样（Lomb-Scargle）
    [psd, f] = plomb(v, t, 'normalized');  % 'normalized' 使单位近似功率/Hz
    title_str = 'PSD (Uneven sampling, Lomb-Scargle, DC excluded)';
else
    % 均匀采样：添加窗函数并用 FFT
    window = hann(N);
    v = v .* window;
    
    yf = fft(v);
    fs = 1 / dt;
    f = (0:floor(N/2)) * (fs / N);
    
    psd = (1 / (fs * N)) * abs(yf(1:floor(N/2)+1)).^2;
    psd(2:end-1) = 2 * psd(2:end-1);
    
    window_correction = N / sum(window.^2);
    psd = psd * window_correction;
    
    title_str = 'PSD (Even sampling, FFT with Hann window, DC excluded)';
end

% 排除 DC 分量（低于0.1 Hz 的近似 DC）
dc_index = find(f < 0.1, 1, 'last');
if ~isempty(dc_index)
    f = f(dc_index+1:end);
    psd = psd(dc_index+1:end);
end

% 如果排除后无数据，停止
if isempty(f)
    error('%s: 排除DC后无数据', csv_path);
end

% 计算 Gamma 频段 (30-100 Hz) 的归一化功率
gamma_mask = (f >= 8) & (f <= 13);
if any(gamma_mask)
    power_gamma = trapz(f(gamma_mask), psd(gamma_mask));  % 使用 trapz 积分功率（更精确）
    power_total = trapz(f, psd);  % 整个频谱总功率 (已排除 DC)
    norm_power_gamma = power_gamma / power_total;  % 归一化功率
else
    norm_power_gamma = 0;  % 如果无 Gamma 频点，设为 0
end

% 计算主频 (PSD 最大值对应的频率)
[~, max_idx] = max(psd);
dominant_freq = f(max_idx);

% 打印结果
fprintf('Gamma 频段归一化功率: %.4f\n', norm_power_gamma);
fprintf('主频: %.2f Hz\n', dominant_freq);

% 绘制频谱图，只显示 0-200 Hz
fig = figure('Visible', 'on');  % 'on' 以显示图表，或 'off' 以不显示
semilogy(f, psd);  % 使用对数 y 轴以提高可读性
xlim([min(f), 200]);
ylim([1e-10, max(psd)*1.1]);
xlabel('Frequency (Hz)');
ylabel('PSD (Power/Hz)');
title(title_str);


fprintf('分析完成！\n');
