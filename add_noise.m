function added_signal = add_noise(signal, snrr, type)
switch (type)
    case 'white'
        added_signal = awgn(signal,snrr,"measured");

    case 'pulse'
        p = 0.0001; %脉冲噪声的概率
        Amax = max(abs(signal))*0.5; %脉冲噪声的最大幅值

            % 生成脉冲噪声掩码（1 表示该点为脉冲噪声）
            pulse_mask = rand(size(signal)) < p;
            % 生成脉冲噪声幅度
            pulse_amplitude = (2 * rand(size(signal)) - 1) .* Amax;
            % 形成稀疏噪声
            pulse_noise = pulse_mask .* pulse_amplitude;
            % 加入信号
            added_signal = signal + pulse_noise;

end
end