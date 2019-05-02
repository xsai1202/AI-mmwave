close all; clear all; clc
%load('BS_codebook_6_6.mat')
%load('MT_codebook_6_6.mat')
load('BS_codebook_4_2.mat');
BS_code_book = code_book;
load('MT_codebook_4_2.mat')
MT_code_book = code_book;
sys.f = 60e9;
sys.lambda = 3e8/sys.f;
sys.phase_hardwire = [0,pi/2,pi,3*pi/2];
si_num =1000;
M=4; N=2;
%% set a phsed array
data = zeros(si_num,64);
data_2 = zeros(si_num,64);
rec_power = zeros(si_num,64);
Sinr = zeros(si_num,64);
rxx = zeros(si_num,2);
obj_location = zeros(si_num,22);
max_ind_curr = zeros(si_num,2);
lineofsight = zeros(si_num,64);
%sys_param = struct('mirrorPlane',[16,50,50,16,50,0,16,-50,0,16,-50,50],'rpos',[],'src',[0,0,8],'obj',[]);
max_ind = zeros(si_num,2);

beam_weight1 = reshape(BS_code_book.beam_weight,[],2,8,1);
beam_weight2 = reshape(MT_code_book.beam_weight,[],2,8,1);

beam_temp_1 = reshape(beam_weight1,8,8);
beam_temp_2 = reshape(beam_weight2,8,8);


for xx = 1:2  
    %tic
    [ rx_1,obj_1 ] = generate_environment_2(3,3,3,3 );
    rx_n = rx_1;
    obj_n = obj_1; 
    if rx_1(2) < -50 
        disp('start')
        for ss = 1:si_num      
            sys_param = struct('mirrorPlane',[16,50,50,16,50,0,16,-50,0,16,-50,50],'rpos',[],'src',[0,0,5],'obj',[]);

            [rx_n, obj_n] = mobilityModule(rx_n,10,0.01,obj_n,3);
            sys_param.rpos = rx_n;
            rxx(ss,:) = rx_n(1:2);  
            obj = cat(1,[-10,0,0,18,60,50],obj_n);
            sys_param.obj = obj;
            obj_location(ss,:) = reshape(obj(2:end,1:2)',1,[]);

            %% set a phased array
            tx_loc = sys_param.src;
            rx_loc = sys_param.rpos;


            %BS_code_book = BS_codebook_6_6.code_book;
            %MT_code_book = MT_codebook_6_6.code_book;
            %% derive the codebook of the phased array
            %beam_weight1 = reshape(BS_code_book.beam_weight,[],2,8,1);
            %beam_weight2 = reshape(MT_code_book.beam_weight,[],2,8,1);
            power_matrix = zeros(8,8);
            power_lineofsight_matrix = zeros(8,8);
            power_db = zeros(8,8);
            power_lineofsight_db = zeros(8,8);
            SINR = zeros(8,8);
            rx_power = zeros(8,8);
            for i =1:8
                for j = 1:8
                    weight_tx = reshape(beam_temp_1(i,:),[],2);
                    weight_rx = reshape(beam_temp_2(j,:),[],2);
                    path_info = getResult(sys,sys_param, weight_rx, weight_tx,M,N);
                    power = calculate_power(path_info,sys);
                    power_lineofsight = calculate_power_lineofsight(path_info,sys,rx_n);
                    rx_power(i,j) = power;
                    noise_p = calculate_noise(path_info); 
                    SINR(i,j) = pow2db(power/noise_p);
                    power_db(i,j) = 10*log10(power)+130;
                    power_lineofsight_db(i,j) = 10*log10(power_lineofsight)+110;
                    power_matrix(i,j) = power_quantization(power_db(i,j),0.1,-100,0);
                    power_lineofsight_matrix(i,j) = power_quantization(power_lineofsight_db(i,j),0.1,-100,0);
                end
            end
            if path_info(1,1).pathgain(1) == 0
                disp(['there exists blockage for situation',num2str(ss)])
            end
            lineofsight(ss,:) = reshape(power_lineofsight_matrix',1,[]);
            data(ss,:) = reshape(power_matrix',1,[]);
            data_2(ss,:) = reshape(power_db',1,[]);
            rec_power(ss,:) = reshape(rx_power',1,[]);
            path_information(ss,:) = reshape(path_info',1,[]);
            Sinr(ss,:) = reshape(SINR',1,[]);
            [a,b] = max(data(ss,:));
%             if n==1
%                 max_ind_curr(ss,:) = [a,b];
%             else
%                 max_ind(ss,:) = [a,b];
%             end
%             disp(['situation ',num2str(ss),' finished'])
            %toc
        end
    end
end
generate_rx_power(rec_power);
generate_SNR(Sinr)
save('data_snrtest.mat')