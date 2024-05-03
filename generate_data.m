function data = generate_data(n,T)

data = generate_exponentially_distributed_data(n,T); % for bi-exponentially distributed data

end

%% function generating bi-exponentially distributed data
function data = generate_exponentially_distributed_data(n,T)
data = zeros(n,T) ;

% generate parameters (a_i's) and lambda

% a      = rand(n,1) ; % the probability of left side of rigin
a = (1/n)*ones(n,1) ;

lambda = rand(n,1) ; % exponential distribution parameter

% generating T samples
for t = 1:T
    data(:,t) = generate_exponential_samples(n,a,lambda); % bi-exponential samples
end
end


%% for bi-exponentially distributed rv's
function sample = generate_exponential_samples(n,a,lambda)
sample = zeros(n,1);

for j = 1:n
    U = rand(1);
    if U < (1 - a(j))
       sample(j) = (1/lambda(j))*log(U/(1 - a(j)));
    else
       sample(j) = (1/lambda(j))*log(a(j)/(1 - U));
    end
end
end
