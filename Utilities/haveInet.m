function tf = haveInet()
%HAVEINET returns true if the owner has internet or false if it hasent
% Taken from
% https://stackoverflow.com/questions/19557118/internet-connection-status-using-matlab
  tf = false;
  try
    java.net.InetAddress.getByName('www.google.com');
    tf = true;
  end
end