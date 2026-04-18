function z = MOP2(x)

    global PROB_NAME;
    prob_name = PROB_NAME;

    if size(x, 1) > 1
        x = x';
    end

    switch prob_name
        %% ==================== DTLZ ====================
        case 'DTLZ1'
            M = 3; 
            D = length(x);
            g = 100 * (D - M + 1 + sum((x(M:end) - 0.5).^2 - cos(20 * pi * (x(M:end) - 0.5))));
            PopObj = 0.5 * repmat(1+g, 1, M) .* fliplr(cumprod([1, x(1:M-1)])) .* [1, 1 - x(M-1:-1:1)];
            z = PopObj';
            
        case 'DTLZ2'
            M = 3;
            g = sum((x(M:end) - 0.5).^2);
            z1 = cos(pi/2 * x(1)) * cos(pi/2 * x(2)) * (1 + g);
            z2 = cos(pi/2 * x(1)) * sin(pi/2 * x(2)) * (1 + g);
            z3 = sin(pi/2 * x(1)) * (1 + g);
            z = [z1; z2; z3];

        case 'DTLZ3'
            M = 3; 
            D = length(x);
            g = 100 * (D - M + 1 + sum((x(M:end) - 0.5).^2 - cos(20 * pi * (x(M:end) - 0.5))));
            z1 = cos(pi/2 * x(1)) * cos(pi/2 * x(2)) * (1 + g);
            z2 = cos(pi/2 * x(1)) * sin(pi/2 * x(2)) * (1 + g);
            z3 = sin(pi/2 * x(1)) * (1 + g);
            z = [z1; z2; z3];
            
        case 'DTLZ4'
            M = 3; 
            alpha = 100;
            g = sum((x(M:end) - 0.5).^2);
            z1 = cos(pi/2 * (x(1)^alpha)) * cos(pi/2 * (x(2)^alpha)) * (1 + g);
            z2 = cos(pi/2 * (x(1)^alpha)) * sin(pi/2 * (x(2)^alpha)) * (1 + g);
            z3 = sin(pi/2 * (x(1)^alpha)) * (1 + g);
            z = [z1; z2; z3];
            
        case 'DTLZ5'
            M = 3;
            g = sum((x(M:end) - 0.5).^2);
            b = x(1);
            c = (1 / (2 + 2*g)) * (1 + 2*g*x(2));
            z1 = cos(pi/2 * b) * cos(pi/2 * c) * (1 + g);
            z2 = cos(pi/2 * b) * sin(pi/2 * c) * (1 + g);
            z3 = sin(pi/2 * b) * (1 + g);
            z = [z1; z2; z3];
            
        case 'DTLZ6'
            M = 3;
            g = sum(x(M:end).^2);
            b = x(1);
            c = (1 / (2 + 2*g)) * (1 + 2*g*x(2));
            z1 = cos(pi/2 * b) * cos(pi/2 * c) * (1 + g);
            z2 = cos(pi/2 * b) * sin(pi/2 * c) * (1 + g);
            z3 = sin(pi/2 * b) * (1 + g);
            z = [z1; z2; z3];
            
        case 'DTLZ7'
            M = 3;
            k = length(x) - M + 1;
            sum_val = sum(x(M:end));
            g = 1 + (9 / k) * sum_val;
            
            z1 = x(1);
            z2 = x(2);
            c1 = (z1 / (1 + g)) * (1 + sin(3 * pi * z1));
            c2 = (z2 / (1 + g)) * (1 + sin(3 * pi * z2));
            h = 3 - (c1 + c2);
            z3 = h * (1 + g);
            z = [z1; z2; z3];

        %% ==================== WFG  ====================
        case 'WFG1'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = ones(1, M-1);
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = zeros(N, K+L);
            t1(:,1:K) = z01(:,1:K); t1(:,K+1:end) = s_linear(z01(:,K+1:end), 0.35);
            t2 = zeros(N, K+L);
            t2(:,1:K) = t1(:,1:K); t2(:,K+1:end) = b_flat(t1(:,K+1:end), 0.8, 0.75, 0.85);
            t3 = b_poly(t2, 0.02);
            t4 = zeros(N, M);
            for i = 1 : M-1
                t4(:,i) = r_sum(t3(:,(i-1)*K/(M-1)+1:i*K/(M-1)), 2*((i-1)*K/(M-1)+1):2:2*i*K/(M-1));
            end
            t4(:,M) = r_sum(t3(:,K+1:K+L), 2*(K+1):2:2*(K+L));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t4(:,M), A(i)) .* (t4(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t4(:,M);
            h = convex(x_wfg); h(:,M) = mixed(x_wfg);    
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        case 'WFG2'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = ones(1, M-1);
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = zeros(N, K+L);
            t1(:,1:K) = z01(:,1:K); t1(:,K+1:end) = s_linear(z01(:,K+1:end), 0.35);
            t2 = zeros(N, K+L/2);
            t2(:,1:K) = t1(:,1:K);
            t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
            t3 = zeros(N, M);
            for i = 1 : M-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)), ones(1,K/(M-1)));
            end
            t3(:,M) = r_sum(t2(:,K+1:K+L/2), ones(1,L/2));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t3(:,M), A(:,i)) .* (t3(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t3(:,M);
            h = convex(x_wfg); h(:,M) = disc(x_wfg);
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        case 'WFG3'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = [1, zeros(1,M-2)];
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = zeros(N, K+L);
            t1(:,1:K) = z01(:,1:K); t1(:,K+1:end) = s_linear(z01(:,K+1:end), 0.35);
            t2 = zeros(N, K+L/2);
            t2(:,1:K) = t1(:,1:K);
            t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
            t3 = zeros(N, M);
            for i = 1 : M-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)), ones(1,K/(M-1)));
            end
            t3(:,M) = r_sum(t2(:,K+1:K+L/2), ones(1,L/2));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t3(:,M), A(:,i)) .* (t3(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t3(:,M);
            h = linear(x_wfg);
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        case 'WFG4'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = ones(1, M-1);
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = s_multi(z01, 30, 10, 0.35);
            t2 = zeros(N, M);
            for i = 1 : M-1
                t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)), ones(1,K/(M-1)));
            end
            t2(:,M) = r_sum(t1(:,K+1:K+L), ones(1,L));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t2(:,M), A(:,i)) .* (t2(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t2(:,M);
            h = concave(x_wfg);
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        case 'WFG5'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = ones(1, M-1);
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = s_decept(z01, 0.35, 0.001, 0.05);
            t2 = zeros(N, M);
            for i = 1 : M-1
                t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)), ones(1,K/(M-1)));
            end
            t2(:,M) = r_sum(t1(:,K+1:K+L), ones(1,L));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t2(:,M), A(:,i)) .* (t2(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t2(:,M);
            h = concave(x_wfg);
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        case 'WFG6'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = ones(1, M-1);
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = zeros(N, K+L);
            t1(:,1:K) = z01(:,1:K); t1(:,K+1:end) = s_linear(z01(:,K+1:end), 0.35);
            t2 = zeros(N, M);
            for i = 1 : M-1
                t2(:,i) = r_nonsep(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)), K/(M-1));
            end
            SUM = zeros(N, 1);
            for i = K+1 : K+L-1
                for j = i+1 : K+L
                    SUM = SUM + abs(t1(:,i)-t1(:,j));
                end
            end
            t2(:,M) = (sum(t1(:,K+1:end),2) + SUM*2) / ceil(L/2) / (1+2*L-2*ceil(L/2));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t2(:,M), A(:,i)) .* (t2(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t2(:,M);
            h = concave(x_wfg);
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        case 'WFG7'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = ones(1, M-1);
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = zeros(N, K+L);
            Y = (fliplr(cumsum(fliplr(z01),2)) - z01) ./ repmat(K+L-1:-1:0, N, 1);
            t1(:,1:K) = z01(:,1:K).^(0.02 + (50-0.02) * (0.98/49.98 - (1-2*Y(:,1:K)) .* abs(floor(0.5-Y(:,1:K))+0.98/49.98)));
            t1(:,K+1:end) = z01(:,K+1:end);
            t2 = zeros(N, K+L);
            t2(:,1:K) = t1(:,1:K); t2(:,K+1:end) = s_linear(t1(:,K+1:end), 0.35);
            t3 = zeros(N, M);
            for i = 1 : M-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)), ones(1,K/(M-1)));
            end
            t3(:,M) = r_sum(t2(:,K+1:K+L), ones(1,L));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t3(:,M), A(:,i)) .* (t3(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t3(:,M);
            h = concave(x_wfg);
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        case 'WFG8'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = ones(1, M-1);
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = zeros(N, K+L);
            t1(:,1:K) = z01(:,1:K);
            Y = (cumsum(z01,2) - z01) ./ repmat(0:K+L-1, N, 1);
            t1(:,K+1:K+L) = z01(:,K+1:K+L).^(0.02 + (50-0.02) * (0.98/49.98 - (1-2*Y(:,K+1:K+L)) .* abs(floor(0.5-Y(:,K+1:K+L))+0.98/49.98))); 
            t2 = zeros(N, K+L);
            t2(:,1:K) = t1(:,1:K); t2(:,K+1:end) = s_linear(t1(:,K+1:end), 0.35);
            t3 = zeros(N, M);
            for i = 1 : M-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)), ones(1,K/(M-1)));
            end
            t3(:,M) = r_sum(t2(:,K+1:K+L), ones(1,L));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t3(:,M), A(:,i)) .* (t3(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t3(:,M);
            h = concave(x_wfg);
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        case 'WFG9'
            PopDec = x; [N,D] = size(PopDec);
            M = 3; K = 2; L = D - K; D_obj = 1; S = 2:2:2*M; A = ones(1, M-1);
            z01 = PopDec ./ repmat(2:2:size(PopDec,2)*2, N, 1);
            t1 = zeros(N, K+L);
            Y = (fliplr(cumsum(fliplr(z01),2)) - z01) ./ repmat(K+L-1:-1:0, N, 1);
            t1(:,1:K+L-1) = z01(:,1:K+L-1).^(0.02 + (50-0.02) * (0.98/49.98 - (1-2*Y(:,1:K+L-1)) .* abs(floor(0.5-Y(:,1:K+L-1))+0.98/49.98)));
            t1(:,end) = z01(:,end);
            t2 = zeros(N, K+L);
            t2(:,1:K) = s_decept(t1(:,1:K), 0.35, 0.001, 0.05);
            t2(:,K+1:end) = s_multi(t1(:,K+1:end), 30, 95, 0.35);
            t3 = zeros(N, M);
            for i = 1 : M-1
                t3(:,i) = r_nonsep(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)), K/(M-1));
            end
            SUM = zeros(N, 1);
            for i = K+1 : K+L-1
                for j = i+1 : K+L
                    SUM = SUM + abs(t2(:,i)-t2(:,j));
                end
            end
            t3(:,M) = (sum(t2(:,K+1:end),2) + SUM*2) / ceil(L/2) / (1+2*L-2*ceil(L/2));
            x_wfg = zeros(N, M);
            for i = 1 : M-1
                x_wfg(:,i) = max(t3(:,M), A(:,i)) .* (t3(:,i)-0.5) + 0.5;
            end
            x_wfg(:,M) = t3(:,M);
            h = concave(x_wfg);
            PopObj = repmat(D_obj*x_wfg(:,M),1,M) + repmat(S,N,1).*h;
            z = PopObj';

        %% ==================== Viennet ====================
        case 'Viennet1'
            z1 = (x(1))^2 + (x(2)-1)^2;
            z2 = (x(1))^2 + ((x(2)+1)^2) + 1;
            z3 = (x(1)-1)^2 + ((x(2))^2) + 2;
            z = [z1; z2; z3];

        case 'Viennet2'
            z1 = 0.5 * (x(1)-2)^2 + 1/13 * (x(2)+1)^2 + 3;
            z2 = (1/36) * ((x(1)+x(2)-3)^2) + ((1/8) * (x(2)-x(1)+2)^2) - 17;
            z3 = (1/175) * ((x(1)+2*x(2)-1)^2) + ((1/17) * (2*x(2)-x(1))^2) - 13;
            z = [z1; z2; z3];

        case 'Viennet3'
            z1 = 0.5 * ((x(1))^2 + (x(2))^2) + sin((x(1))^2 + (x(2))^2);
            z2 = ((3*x(1)-2*x(2)+4)^2)/8 + ((x(1)-x(2)+1)^2)/27 + 15;
            z3 = (1/((x(1))^2+(x(2))^2+1)) - 1.1 * exp(-(x(1))^2 - (x(2))^2);
            z = [z1; z2; z3];

        otherwise
            error([ prob_name, ]);
    end
end

%% ========================================================================
% ========================================================================
function Output = s_linear(y,A)
    Output = abs(y-A) ./ abs(floor(A-y)+A);
end

function Output = b_flat(y,A,B,C)
    Output = A + min(0,floor(y-B))*A.*(B-y)/B - min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
    Output = round(Output*1e4)/1e4;
end

function Output = b_poly(y,a)
    Output = y.^a;
end

function Output = r_sum(y,w)
    Output = sum(y .* repmat(w, size(y,1), 1), 2) ./ sum(w);
end

function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1), 1-cos(x(:,1:end-1)*pi/2)], 2)) .* [ones(size(x,1),1), 1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
    Output = 1 - x(:,1) - cos(10*pi*x(:,1)+pi/2)/10/pi;
end

function Output = disc(x)
    Output = 1 - x(:,1) .* (cos(5*pi*x(:,1))).^2;
end

function Output = linear(x)
    Output = fliplr(cumprod([ones(size(x,1),1), x(:,1:end-1)], 2)) .* [ones(size(x,1),1), 1-x(:,end-1:-1:1)];
end

function Output = s_multi(y,A,B,C)
    Output = (1 + cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C))) + 4*B*(abs(y-C)/2./(floor(C-y)+C)).^2) / (B+2);
end

function Output = concave(x)
    Output = fliplr(cumprod([ones(size(x,1),1), sin(x(:,1:end-1)*pi/2)], 2)) .* [ones(size(x,1),1), cos(x(:,end-1:-1:1)*pi/2)];
end

function Output = s_decept(y,A,B,C)
    Output = 1 + (abs(y-A)-B) .* (floor(y-A+B)*(1-C+(A-B)/B)/(A-B) + floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B) + 1/B);
end

function Output = r_nonsep(y,A)
    Output = zeros(size(y,1), 1);
    for j = 1 : size(y,2)
        Temp = zeros(size(y,1), 1);
        for k = 0 : A-2
            Temp = Temp + abs(y(:,j) - y(:,1+mod(j+k,size(y,2))));
        end
        Output = Output + y(:,j) + Temp;
    end
    Output = Output ./ (size(y,2)/A) / ceil(A/2) / (1+2*A-2*ceil(A/2));
end