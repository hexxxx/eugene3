Rails.application.routes.draw do
  root to: 'proteins#index'
  resources :proteins
end
