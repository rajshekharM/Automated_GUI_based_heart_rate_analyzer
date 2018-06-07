function video_cropping(video_file,start_frame,stop_frame,L_user)

start_frame_user = str2double(start_frame);
stop_frame_user = str2double(stop_frame);
L = str2double(L_user);

[filepath,name,ext] = fileparts(video_file)
input_video=strcat(name,ext)
transcoded_segment_of_input_video='time_cropped_video.avi';
fps=30.0  % desired constant frame per second rate (CFR)
cd (filepath)
cd
obj=VideoReader(video_file);
start_frame = max (start_frame_user - (L-1)/2, 1)
stop_frame = min (stop_frame_user + (L-1)/2, obj.NumberOfFrames )
stop_frame=max(stop_frame,start_frame);

%put the file path(after 'cd ') where HnadbrakeeCLI.exe is located \\ the
%videos should also be located/copied to the same path
s=sprintf('cd %s',filepath);
system(s);
    
s=sprintf('HandBrakeCLI.exe -i %s --start-at frame:%d  --stop-at frame:%d -e x264 --encoder-preset veryfast -q 10 --cfr -r %d -o %s', input_video, start_frame, stop_frame-start_frame+1,fps,transcoded_segment_of_input_video);
system(s);

obj=VideoReader(transcoded_segment_of_input_video);
for k=1:obj.NumberOfFrames
  mov(k).cdata = read(obj, k);
end
%write each frame from mov array into video with start stop frames
%specified

%point to the start frame of transcoded(time-cropped)video for image cropping tool to work
k=1;
I=mov(k).cdata;

vidObj_crop = VideoWriter('final_video.avi');
open(vidObj_crop);
% Give the user to make the rectangle on the image to Crop image starting 
[J, rect] = imcrop(I);

%Apply the same rectangle to rest of the frames of the whole video 
for k=1:obj.NumberOfFrames-1   %it reads one frame more
  I = read(obj, k);
  cropped_img = imcrop(I, rect);
  writeVideo(vidObj_crop,cropped_img);
  display(['frames cropped: ' num2str(k+start_frame)]);
end

close(vidObj_crop);