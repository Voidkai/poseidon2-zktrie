use lazy_static::lazy_static;
use std::sync::Arc;

use crate::{
    fields::bn256::FpBN256,
    neptune::neptune_params::NeptuneParams,
};

lazy_static! {
    // Number of partial rounds:
    // ceil(1.125 * ceil(log_d(2) * (min(kappa, log_2(p)) - 6) + 3 + t + log_d(t)))
    // BN256
    pub static ref NEPTUNE_BN_PARAMS: Arc<NeptuneParams<FpBN256>> = Arc::new(NeptuneParams::new(4, 5, 6, 68));
}