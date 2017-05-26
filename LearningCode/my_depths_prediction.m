

function my_depths_prediction()

warning off
run( 'InitialPath(false)');


opts=[];

% img_dir='../../Dataset/Make3D/';
% img_type='outdoor';
% out_dir = '../../result/Make3D/';

img_dir='../../Dataset/NYUD2/';
img_type='indoor';
out_dir = '../../result/NYUD2/';


% if you want to try your own images,
% create a dataset name, put your images in '<img_dir>/',
% then specify the trained model.
% choose indoor or outdoor scene trained model, depends on your images.

% % some indoor image examples:
% img_dir='../custom_indoor_sample/'; 
% img_type='indoor';


% some outdoor image examples:
% img_dir='../custom_outdoor_sample/';
% img_type='outdoor';



%folder for saving prediction results
[tmp_dir, img_dir_name]=fileparts(img_dir);
if isempty(img_dir_name)
    [~, img_dir_name]=fileparts(tmp_dir);
end
result_dir=fullfile('..', '..', 'results', img_dir_name);


file_infos=dir(fullfile(img_dir,'*'));
valid_file_flags=true(length(file_infos), 1);
for f_idx=1:length(file_infos)
    if file_infos(f_idx).isdir
        valid_file_flags(f_idx)=false;
    end
end
file_infos=file_infos(valid_file_flags);

img_num=length(file_infos);

if img_num==0
    error('Error! no test images found!');
end

err = zeros(img_num, 3);
accu = zeros(img_num, 3);
for img_idx=1:img_num
    
    one_f_info=file_infos(img_idx);
    one_img_file=one_f_info.name;
    full_img_file=fullfile(img_dir, one_img_file);

    fprintf('\n-------------------------------------------\n');
    fprintf('processing image (%d of %d): %s\n', img_idx, img_num, full_img_file);
 
    
    if strcmpi(img_type, 'outdoor')
        label_dir = '../../Dataset/Make3D_data/Gridlaserdata/';
        temp_id = strfind(one_img_file, '.jpg');
        label_file = strcat('depth_sph_corr', one_img_file(4:temp_id), 'mat');
        full_label_file = fullfile(label_dir, label_file);
    end

    if strcmpi(img_type, 'indoor')
        label_dir = '../../Dataset/NYUD2_data/nyu_test_label/';
        temp_id = strfind(one_img_file, '.png');
        label_file = strcat(one_img_file(1:temp_id), 'mat');
        full_label_file = fullfile(label_dir, label_file);
    end

    output_dir = strcat(out_dir, one_img_file(1:temp_id-1), '/');
    depth_pre = OneShot3dEfficient(full_img_file, output_dir);    
    
    load(full_label_file);
    if strcmpi(img_type, 'outdoor')
        ground_truth = Position3DGrid(:, :, 4);
    end

    if strcmpi(img_type, 'indoor')
        ground_truth = depths;
    end
    
    thr = 1.25;
    accu(img_idx, :) = do_accuracy_evaluate(depth_pre, ground_truth, thr);
    err(img_idx, :) = do_error_evaluate(depth_pre, ground_truth);
    close all
    
end
err_ave = sum(err) / img_num;
accu_ave = sum(accu) / img_num;
fprintf('%s%d\n%s%d\n%s%d\n', 'rel:', err_ave(1), 'rms:', err_ave(2), ...
    'log10:', err_ave(3));
fprintf('%s\n', '------------------------------------');
fprintf('%s%d\n%s%d\n%s%d\n', 'accu_1.25:', accu_ave(1), 'accu_1.25^2:', ...
    accu_ave(2), 'accu_1.25^3:', accu_ave(3));


end

